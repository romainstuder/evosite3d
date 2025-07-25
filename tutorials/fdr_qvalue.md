False discovery rate correction for multiple testing of positive selection [Pratictal]

When performing a test for the detection of positive selection (i.e. with CodeML/PAML) on multiple
branches one by one, we need to correct your result to avoid false discovery rate (there are many
methods (see the wiki page: <http://en.wikipedia.org/wiki/False_discovery_rate>).

One I like is the QVALUE package, for R:
<http://www.bioconductor.org/packages/release/bioc/html/qvalue.html>

Launch R and install QVALUE by typing:

```shell
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
```

` and load it:

```shell
library(qvalue)
```

We need to load the list of p-value into R. This is a single with only p-values. The order is really
important, so we have in a separate file the ordered names of the branch and/or gene tested.

```shell

0.00623359809654
0.148525646364
1.0
0.0682850741607
0.0322125346865
```

A file containing p-values from 760 gene families is provided here as an example:
<https://docs.google.com/file/d/0BxcXpZeUylGbMFd1emJ1MERKems/edit?usp=sharing>

Once we have our file of p-values, we can load it into R:

```shell
p <- scan("pvalues.list", na.strings=T)
```

First, have a look at the distribution:

hist(p, breaks = 20, main = paste("Distribution of p-values"), xlab="Value")

which will display this histogram: ￼ The distribution is bimodal, as expected: when the test is
negative, we will have a p-value = 1, and when the test is significant, the p-value is generally
very low.

Now, we can perform the qvalue correction: qobj <- qvalue(p, pi0.meth="bootstrap", fdr.level=0.05)

Two comments:

1. The option “bootstrap” is really important, as the distribution is bimodal (see page 11 of the
   QVALUE manual).
2. The fdr.level=0.05 is a threshold which estimates that 5% of our tests could be false-positive.

Finally, we write the output into a file: write.qvalue(qobj, file="qvalues.list")

The output file has five columns:

1. The p-values, in the same order as in the input file.
2. The corresponding q-value.
3. The pi0 value.
4. The significance of the test, based on the threshold (0 or 1).
5. The fdr threshold used to assign as significant.

```
"pvalue" "qvalue" "lfdr" "pi0" "significant" "fdr.level"
0.00623359809654 0.00598773686818557 0.0291702400633863 0.323886639676113 1 0.05
0.148525646364 0.0731203182099692 0.391218512376612 0.323886639676113 0 0.05
1 0.323886639676113 1 0.323886639676113 0 0.05
0.0682850741607 0.0386405371024297 0.206113388985345 0.323886639676113 1 0.05
0.0322125346865 0.0212580142290782 0.111050362300586 0.323886639676113 1 0.05
```

Now we can say which test is significant and which is not.

Please don't hesitate to contact me if you have any comments or questions.
