version 1.0

task quality_check {
    input {
        File data
        String prefix
        String docker = "docker.io/library/singlecell:1.0.0"
        Int memory_gb = 16
        Int cpu = 16
    }
    Int disk_size = ceil(size(data, "GiB")) + 2
    command <<<
        ln -s ~{data} ~{prefix}.tar.gz
        tar -xzvf ~{prefix}.tar.gz && rm ~{prefix}.tar.gz
        Rscript /scripts/qc_detect_doublets.R -s ~{prefix}
        tar -czvf ~{prefix}_qc.tar.gz *_passed_QC_noDoublet.txt 
    >>>
    output {
        File out = prefix + "_qc.tar.gz"
    }
    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
    }
}