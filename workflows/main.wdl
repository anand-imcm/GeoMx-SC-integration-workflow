version 1.0

import "./tasks/extract_data.wdl" as ex
import "./tasks/quality_check.wdl" as qc
import "./tasks/enrichment.wdl" as en

workflow main {
    input {
        File dataset
        String prefix
        File geneset
        File metadata
    }
    call ex.extract_data {
        input: data = dataset, prefix = prefix
    }
    scatter (sample_tar in extract_data.samples_tar) {
        call qc.quality_check {
            input: data = sample_tar
        }
    }
    call en.enrichment {
        input: data = quality_check.out, prefix = prefix, geneset = geneset, metadata = metadata
    }
    output {
        Array[File] qc_matrix = quality_check.out
    }
}