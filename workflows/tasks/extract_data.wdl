version 1.0

task extract_data {
    input {
        File data
        String prefix
        String docker = "ubuntu:20.04"
        Int memory_gb = 16
        Int cpu = 16
    }
    Int disk_size = ceil(size(data, "GiB")* 10) + 10
    command <<<
        ln -s ~{data} ~{prefix}.tar
        tar -xvf ~{prefix}.tar && rm ~{prefix}.tar
        find . -type f -name "*.tar.gz" -exec mv {} $(pwd) \;
    >>>
    output {
        Array[File] samples_tar = glob("*.tar.gz")
    }
    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
    }
}