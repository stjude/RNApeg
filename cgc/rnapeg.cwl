{
    "class": "CommandLineTool",
    "cwlVersion": "v1.1",
    "baseCommand": [],
    "inputs": [
        {
            "id": "fasta",
            "type": "File",
            "inputBinding": {
                "prefix": "-f",
                "shellQuote": false,
                "position": 32
            },
            "label": "FASTA File",
            "doc": "FASTA file",
            "sbg:fileTypes": "FA",
            "secondaryFiles": [
                {
                    "pattern": ".fai"
                }
            ]
        },
        {
            "id": "bam_file",
            "type": "File",
            "inputBinding": {
                "prefix": "-b",
                "shellQuote": false,
                "position": 2
            },
            "label": "BAM File",
            "doc": "input bamfile mapped to human genome builds GRCh37-lite or GRCh38_no_alt.",
            "sbg:fileTypes": "BAM",
            "secondaryFiles": [
                {
                    "pattern": ".bai?"
                },
                {
                    "pattern": "$(self.nameroot).bai?"
                }
            ]
        },
        {
            "id": "refFlat",
            "type": "File",
            "inputBinding": {
                "prefix": "-r",
                "shellQuote": false,
                "position": 50
            },
            "label": "refFlat Text",
            "doc": "Uncompressed refFlat.txt file from the UCSC genome annotation database",
            "sbg:fileTypes": "TXT"
        }
    ],
    "outputs": [
        {
            "id": "output_file",
            "type": "File",
            "outputBinding": {
                "glob": "$(inputs.bam_file.nameroot).bam.junctions.tab.shifted.tab"
            }
        }
    ],
    "label": "rnapeg",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "ResourceRequirement",
            "ramMin": 4000,
            "coresMin": 0
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "hints": [
        {
            "class": "DockerRequirement",
            "dockerPull": "ghcr.io/stjude/rnapeg:latest"
        }
    ]
}
