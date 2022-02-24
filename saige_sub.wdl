task test {

    File nullfile
    File varianceratiofile = sub(nullfile, ".rda", ".varianceRatio.txt")
    String outfileprefix = basename(nullfile, ".rda") + "-"
    Array[File] bgenfiles
    File samplefile
    Int minmac
    String docker
    String loco
    String analysisType
    String cohort

    command {

        python3 <<EOF
        import os
        import subprocess
        import time
        processes = set()
        cmd_prefix = 'export MKL_NUM_THREADS=1; export MKL_DYNAMIC=false; export OMP_NUM_THREADS=1; export OMP_DYNAMIC=false; \
            step2_SPAtests.R \
            --minMAC=${minmac} \
            --sampleFile=${samplefile} \
            --GMMATmodelFile=${nullfile} \
            --varianceRatioFile=${varianceratiofile} \
            --numLinesOutput=1000 \
            --IsOutputAFinCaseCtrl=TRUE \
            --IsOutputNinCaseCtrl=TRUE \
            --IsOutputHetHomCountsinCaseCtrl=TRUE \
            --LOCO=${loco} '
        for file in '${sep=" " bgenfiles}'.split(' '):
            cmd = cmd_prefix + '--bgenFile=' + file
            if ("X" in file):
                cmd = cmd + ' --is_rewrite_XnonPAR_forMales=TRUE \
                              --X_PARregion=10001-2781479,155701383-156030895'
            cmd = cmd + ' --SAIGEOutputFile=${outfileprefix}' + os.path.basename(file) + '.SAIGE.txt'
            logfile = open('SAIGE_log_${outfileprefix}' + os.path.basename(file) + '.txt', 'w')
            processes.add(subprocess.Popen(cmd, shell=True, stdout=logfile))
        print(time.strftime("%Y/%m/%d %H:%M:%S") + ' ' + str(len(processes)) + ' processes started', flush=True)
        n_rc0 = 0
        while n_rc0 < len(processes):
            time.sleep(60)
            n_rc0 = 0
            for p in processes:
                p_poll = p.poll()
                if p_poll is not None and p_poll > 0:
                    raise('subprocess returned ' + str(p_poll))
                if p_poll == 0:
                    n_rc0 = n_rc0 + 1
            print(time.strftime("%Y/%m/%d %H:%M:%S") + ' ' + str(n_rc0) + ' processes finished', flush=True)
        EOF
    }

	output {
        Array[File] out = glob("*.SAIGE.txt")
        Array[File] logs = glob("SAIGE_log_*.txt")
    }

    runtime {
        docker: "${docker}"
        cpu: length(bgenfiles)
        memory: (4 * length(bgenfiles)) + " GB"
        disks: "local-disk " + (length(bgenfiles) * ceil(size(bgenfiles[0], "G")) + 1) + " HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}


task combine {

    String pheno
    String traitType
    Array[Array[File]] results2D
    Array[File] results = flatten(results2D)
    File blacklist
    Float info_threshold
    String chrcol
    String p_valcol
    String bp_col
    Int loglog_pval
    String docker
    String cohort

    command <<<

        for file in ${sep=" " results}; do
            if [[ $file == *.gz ]]
            then
                gunzip -c $file > $file"DATAUNZIP"
            else
                mv $file `basename $file`"DATAUNZIP"
            fi
        done

        if [[ ${traitType} == "binary" ]]; then
            cat <(head -n 1 `basename ${results[0]}"DATAUNZIP"` | tr ' ' '\t') \
            <(awk 'FNR>1 { printf "%s\t%d\t%s\t%s\t%s\t%s\t%.2f\t%.3e\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.3e\t%.3e\t%d\t%.4f\t%.3e\t%.3e\t%.3e\t%d\t%d\t%.3e\t%.3e\t%.3e\t%.3e\n", \
            $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26 }' \
            `find *DATAUNZIP | sort -V | tr '\n' ' '`) | sort -k1,1V -k2,2g -s > ${pheno}.${cohort}
        else
            cat <(head -n 1 `basename ${results[0]}"DATAUNZIP"` | tr ' ' '\t') \
            <(awk 'FNR>1 { printf "%s\t%d\t%s\t%s\t%s\t%s\t%.2f\t%.3e\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.3e\t%.3e\t%d\t%.3e\t%.3e\n", \
            $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22 }' \
            `find *DATAUNZIP | sort -V | tr '\n' ' '`) | sort -k1,1V -k2,2g -s > ${pheno}.${cohort}
        fi
        postscript.py ${pheno}.${cohort} ${blacklist} ${info_threshold} > ${pheno}.${cohort}.pheweb && \
        qqplot.R --file ${pheno}.${cohort}.pheweb --bp_col "${bp_col}" --pval_col "${p_valcol}" --chrcol "${chrcol}" --loglog_pval ${loglog_pval} && \
        bgzip ${pheno}.${cohort} && \
        tabix -S 1 -b 2 -e 2 -s 1 ${pheno}.${cohort}.gz && \
        bgzip ${pheno}.${cohort}.pheweb && \
        tabix -S 1 -b 2 -e 2 -s 1 ${pheno}.${cohort}.pheweb.gz

    >>>

    output {
        File out = pheno + "." + cohort + ".gz"
        File out_ind = pheno + "." + cohort + ".gz.tbi"
        File out_pheweb = pheno + "." + cohort + ".pheweb.gz"
        File out_pheweb_ind = pheno + "." + cohort + ".pheweb.gz.tbi"
        File qq = pheno + "." + cohort + ".pheweb_" + p_valcol + "_qqplot.png"
        File manh = pheno + "." + cohort + ".pheweb_" + p_valcol + "_manhattan.png"
        File manh_loglog = pheno + "." + cohort + ".pheweb_" + p_valcol + "_manhattan_loglog.png"
        File quantiles = pheno + "." + cohort + ".pheweb_" + p_valcol + "_qquantiles.txt"
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "30 GB"
        disks: "local-disk 100 HDD"
        zones: "europe-west1-b"
        preemptible: 1
        noAddress: true
    }
}

workflow test_combine {

    String docker
    String pheno
    String traitType
    String nullfile
    File bgenlistfile
    Array[Array[String]] bgenfiles2D = read_tsv(bgenlistfile)
    String loco
    String analysisType
    String cohort

    scatter (bgenfiles in bgenfiles2D) {
        call test {
            input: docker=docker, nullfile=nullfile, bgenfiles=bgenfiles, loco=loco, analysisType=analysisType, cohort=cohort
        }
    }

    call combine {
        input: pheno=pheno, traitType=traitType, results2D=test.out, docker=docker, cohort=cohort
    }

    output {
        File out = combine.out
    }
}