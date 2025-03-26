nextflow run nf-core/rnaseq \
  -r 3.18.0 \
  -params-file rnaseq_params.yml \
  -profile docker \
  -resume
