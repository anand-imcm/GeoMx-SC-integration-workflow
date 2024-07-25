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
        File matrix_qc = quality_check.matrix_qc
        File cells_qc = quality_check.cells_qc
    }
}