version 1.0

task quality_check {
    input {
        File data
        String prefix
        String docker = "docker.io/library/singlecell:1.0.0"
        Int memory_gb = 16
        Int cpu = 16
    }
    # Int disk_size_gb = ceil(size(data, "GB")) + 2
    command <<<
        ln -s ~{data} ~{prefix}.tar.gz
        tar -xzvf ~{prefix}.tar.gz
        ls .
        Rscript /scripts/qc_detect_doublets.R -s ~{prefix}
    >>>
    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory_gb}GB"
        # disks: "local-disk ~{disk_size_gb} HDD"
    }
}