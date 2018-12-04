#!/usr/bin/env cwl-runner

class: CommandLineTool

id: "ppcg-qc-from-sanger"

label: "tool to extract PPCG defined QC metrics"

cwlVersion: v1.0

doc: |
    ![build_status](https://quay.io/repository/wtsicgp/ppcg-qc-from-sanger/status)
    A Docker container of ppcg-qc-from-sanger. See the [ppcg-qc-from-sanger](https://github.com/cancerit/ppcg-qc-from-sanger) website for more information.

dct:creator:
  "@id": "yaobo.xu@sanger.ac.uk"
  foaf:name: Yaobo Xu
  foaf:mbox: "yx2@sanger.ac.uk"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/wtsicgp/ppcg-qc-from-sanger:0.2.4"

inputs:
  tumour_bas:
    type:
      type: array
      items: File
      inputBinding:
        prefix: -tb
        shellQuote: true
        separate: true
    doc: "tumour sample BAS file from cgpmap pipeline"

  normal_bas:
    type:
      type: array
      items: File
      inputBinding:
        prefix: -nb
        shellQuote: true
        separate: true
    doc: "normal sample BAS file from cgpmap pipeline"

  metadata:
    type:
      - "null"
      - type: array
        items: File
        inputBinding:
          prefix: -mt
          shellQuote: true
          separate: true
    doc: "metadata files"

  cgpwgs_result_tar:
    type:
      type: array
      items: File
      inputBinding:
        prefix: -rt
        shellQuote: true
        separate: true
    doc: "variant calling results tar file from cgpwgs"

  output_tar_name:
    type: string
    default: out.tar.gz
    doc: "output tar file name"
    inputBinding:
      prefix: -o
      position: 4
      separate: true

  count_variants:
    type: boolean?
    doc: "output variant counts"
    inputBinding:
      prefix: -cv

outputs:
  result_tar:
    type: File
    outputBinding:
      glob: $(inputs.output_tar_name)

baseCommand: ["ppcg-qc-from-sanger"]
