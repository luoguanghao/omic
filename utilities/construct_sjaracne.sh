
# example
sudo sjaracne local \
-e ../inputs/BRCA100.exp \
-g ../inputs/tf.txt  \
-n 2 \
-o ./SJARACNe_output.final


# NetBID2 Nat. Cancer
sjaracne local -e ./nodup_ALL_T_264_sample_expr_mat.txt \
-g ./tf_symbol.txt \
-n 2 \
-o ./SJARACNE_out.final

# VIZOME



pip install numpy==1.17.0
pip install scipy==1.0.1
pip install pandas==0.22.0



java -Xmx5G -jar ../dist/aracne.jar -e ./matrix.txt  -o ./ --tfs ./tfs.txt --pvalue 1E-8 --seed 1 \
--calculateThreshold



java -Xmx5G -jar ../dist/aracne.jar -e ./matrix.txt  -o ./ --tfs ./tfs.txt --pvalue 1E-8 --seed 1 \
--nobootstrap --noDPI

java -Xmx5G -jar ../dist/aracne.jar -e ./matrix.txt  -o ./ --tfs ./tfs.txt --pvalue 1E-8 --seed 1




