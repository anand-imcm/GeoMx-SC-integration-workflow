version 1.0

task quality_check {
    input {
        File data
        String docker = "docker.io/library/singlecell:1.0.0"
        Int memory_gb = 16
        Int cpu = 16
    }
    Int disk_size = ceil(size(data, "GiB")) + 2
    command <<<
        filename=$(basename ~{data})
        base="${filename%%_*}" 
        ln -s ~{data} ${base}.tar.gz
        tar -xzvf ${base}.tar.gz && rm ${base}.tar.gz
        Rscript /scripts/qc_detect_doublets.R -s ${base}
        tar -czvf ${base}.tar.gz *_passed_QC_noDoublet.txt
        echo "${base}.tar.gz" > output_filename.txt
    >>>
    output {
        File out = read_string("output_filename.txt")
    }
    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
    }
}