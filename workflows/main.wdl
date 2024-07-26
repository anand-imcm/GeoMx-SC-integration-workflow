version 1.0

import "./tasks/extract_data.wdl" as ex
import "./tasks/quality_check.wdl" as qc

workflow main {
    input {
        File dataset
        String prefix
    }
    call ex.extract_data {
        input: data = dataset, prefix = prefix
    }
    scatter (sample_tar in extract_data.samples_tar) {
        call qc.quality_check {
            input: data = sample_tar
        }
    }

    output {
        Array[File] qc_matrix = quality_check.out
    }
}