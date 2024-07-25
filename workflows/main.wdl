version 1.0

import "./tasks/quality_check.wdl" as qc

workflow main {
    input {
        File dataset
        String prefix
    }
    call qc.quality_check {
        input: data = dataset, prefix = prefix
    }
    output {
        File qc_matrix = quality_check.out
    }
}