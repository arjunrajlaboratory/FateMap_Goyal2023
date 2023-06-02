###
Alterations to code:

2019/08/06
- running on MAC OSX "zcat" doesn't work. You need to use "zcat <" or "gunzip -c".
  this makes it difficult to directly implement the wrapper. to deal with this I changed the
  rules/align.cmk and rules/trim.cmk to have input and output that aren't gzipped.
