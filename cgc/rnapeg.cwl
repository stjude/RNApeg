{
    "class": "CommandLineTool",
    "cwlVersion": "v1.2",
    "$namespaces": {
        "sbg": "https://sevenbridges.com"
    },
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
            "doc": "FASTA file of reference sequences. Must match the reference to which the BAM is aligned.",
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
            "doc": "Input BAM file mapped to human genome builds GRCh37-lite or GRCh38_no_alt.",
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
    "doc": "RNApeg is an RNA junction calling, correction, and quality-control package.\n\n## Inputs\n* **Fasta** - Reference genome in FASTA format. Chromosomes must match those in BAM header.\n* **BAM** - Aligned RNA-Seq BAM\n* **Refflat** - Uncompressed annotation file from UCSC genome annotation database\n\n## Outputs\n* **Junctions** - Flat file containing called RNA junctions",
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
            "dockerPull": "cgc-images.sbgenomics.com/stjude/rnapeg:latest"
        }
    ],
    "sbg:links": [
        {
            "id": "https://rnajournal.cshlp.org/content/24/8/1056.short",
            "label": "Publication"
        },
        {
            "id": "https://github.com/stjude/RNApeg",
            "label": "Source Code"
        }
    ],
    "sbg:wrapperLicense": "Apache 2.0 License",
    "sbg:categories": [
        "RNA-Seq",
        "Junction Calling"
    ],
    "sbg:license": "Apache 2.0 License"
}
