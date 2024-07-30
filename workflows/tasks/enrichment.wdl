version 1.0

task enrichment {
    input {
        Array[File] data
        File geneset
        File metadata
        String prefix
        String docker = "docker.io/library/singlecell:1.0.0"
        Int memory_gb = 16
        Int cpu = 16
    }
    Int disk_size = ceil(size(data, "GiB")) + 2
    command <<<
        for file in ~{sep=' ' data}; do
            filename=$(basename ${file})
            base="${filename%%.*}"
            mkdir -p ${base}/QC
            tar -xzvf ${file} -C ${base}/QC
            echo ${base}
        done
        Rscript /scripts/enrichment_analysis.R -o . -s ~{prefix} -g ~{geneset} -m ~{metadata}
    >>>

    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
    }
}