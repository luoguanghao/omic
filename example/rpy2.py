

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

from rpy2.robjects.conversion import localconverter


# import pandas.DataFrame:Dactolisib_CCLE_sen_unsen_expr_mat into r data.frame object
with localconverter(ro.default_converter + pandas2ri.converter):
    r_Dactolisib_CCLE_sen_unsen_expr_mat = ro.conversion.py2rpy(Dactolisib_CCLE_sen_unsen_expr_mat)

# import list/array into r object
r_condition_table = robjects.StrVector(condition_table)

# 执行一段代码
rcode = """
expr_mat <- %s
condition_table <- %s

dds <- DESeqDataSetFromMatrix(expr_mat, DataFrame(condition_table), design= ~ condition_table )
dds2 <- DESeq(dds)
res <- results(dds2)
"""%(r_CCLE_sen_unsen_expr_mat.r_repr(), r_condition_table.r_repr())

robjects.r(rcode)


# python下生成的robject导入robjects.r
robjects.r(r_condition_table.r_repr())

# 直接使用R内部函数
robjects.r('pi')










