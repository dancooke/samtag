# samtag

[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/dancooke/samtag)
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/dancooke/samtag?label=docker)](https://hub.docker.com/r/dancooke/samtag)

`samtag` is a command-line tool for tagging reads in a BAM/CRAM file whose name appear in an input list. It is possible to specify default a tag/flag applied to every read and/or per-read tags/flags. `samtag` can also generate statistics about tags and their values based on regular expressions.

## Usage

The minimal input is a BAM/CRAM file and a TSV specifying read names in the first column, a tag in the second column (optional), and flags in the third column (optional). The format for tags is `TAG:VALUE`, unless a default tag is provided with the `--tag` option without a value, in which case the `TAG` names in the input file can be ommitted (i.e. just specify `VALUE`).

```shell
# Add "ZA:Z:FOO" tag to all reads in test/reads1.tsv,
# in addition to tags/flags specified in test/reads1.tsv.
$ samtag tag --tag ZA:BAR -o tagged1.bam test/test.bam test/reads1.tsv

# Add "ZA" tag to all reads in test/reads2.tsv,
# with the values specified in test/reads2.tsv.
$ samtag tag --tag ZA -o tagged2.bam test/test.bam test/reads2.tsv
```

To compute statistics about tags:

```shell
$ samtag stats --tag ZA tagged1.bam
tag     value   count   fraction
*       *       3859    1
ZA      *       8       0.00207308

$ samtag stats --tag ZA --split --sort tagged1.bam
tag     value   count   fraction
*       *       3859    1
ZA      *       8       0.00207308
ZA      BAR     7       0.00181394
ZA      FOO     1       0.000259134

$ samtag stats --tag ZA:BAR tagged1.bam
tag     value   count   fraction
*       *       3859    1
ZA      BAR     7       0.00181394
```
