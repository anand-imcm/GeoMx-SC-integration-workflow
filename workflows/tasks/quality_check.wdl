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
        tar -xzvf ~{prefix}.tar.gz
        Rscript /scripts/qc_detect_doublets.R -s ~{prefix}
    >>>
    output {
        File matrix_qc = prefix + "_Gene_Expression_Matrix_passed_QC_noDoublet.txt"
        File cells_qc = prefix + "_Cells_passed_QC_noDoublet.txt"
    }
    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
    }
}