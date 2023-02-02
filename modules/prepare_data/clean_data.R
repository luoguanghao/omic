
##################
# 检查数据是否都合法
##################





# 包括各种转换，tpm转rpkm等等
#
#
#

###########
# rpkm转tpm
###########
# raw_exp_mat: gene*sample rpkm without log
## out:
# exp_mat_tpm gene*sample, tpm without log, sum of every columns is 1e6
# ===
fpkmToTpm <- function(fpkm)
{
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
exp_mat_tpm = apply(raw_exp_mat,2,fpkmToTpm)















