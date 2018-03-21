# vest-transform

## setup

### install our dependencies
pip install -r requirements.txt

### data
#### initialize VMC database

See https://github.com/biocommons/biocommons.seqrepo/blob/master/doc/mirror.rst#fetching-using-rsync-manually

I rsync'd the entire repo, set instance name to master

export SEQREPO_INSTANCE_NAME=master
export SEQREPO_ROOT_DIR=/tmp/seqrepo

#### scores
untar swift://biostream/source/vest-3.0-precompute.tar.gz


## output

```
$ head -1  packaged-scores/chr1.all.scores
chr1	35143	T	G	T	P	0.315	ENST00000417324:T85P

$ python transform.py
# raw transform
{
  "aa_mutation": "p.T85P",
  "alt": "G",
  "chromosome": "1",
  "protein_alt": "P",
  "protein_ref": "T",
  "ref": "T",
  "score": 0.315,
  "start": 35143,
  "transcript": "ENST00000417324"
}
# VMC identifiers
{
  "alleles": {
    "VMC:GA_Wu1Zi_4-67-c6qYnYHywtcdEW5AoY2b_": {
      "id": "VMC:GA_Wu1Zi_4-67-c6qYnYHywtcdEW5AoY2b_",
      "location_id": "VMC:GL_ar0lWY-kTvmNuB5FtJK7bCJX_7XuTuwV",
      "state": "G"
    }
  },
  "identifiers": {
    "VMC:GS_S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU": [
      {
        "accession": "NC_000001.10",
        "namespace": "NCBI"
      }
    ]
  },
  "locations": {
    "VMC:GL_ar0lWY-kTvmNuB5FtJK7bCJX_7XuTuwV": {
      "id": "VMC:GL_ar0lWY-kTvmNuB5FtJK7bCJX_7XuTuwV",
      "interval": {
        "end": 35143,
        "start": 35143
      },
      "sequence_id": "VMC:GS_S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU"
    }
  },
  "meta": {
    "aa_mutation": "p.T85P",
    "transcript": "ENST00000417324",
    "vest_score": 0.315,
    "vmc_version": 0
  }
}
```

see https://docs.google.com/document/d/12E8WbQlvfZWk5NrxwLytmympPby6vsv60RxCeD5wc1E/edit#heading=h.ite1lnpy1bca for more
