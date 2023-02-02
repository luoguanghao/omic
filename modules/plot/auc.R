





roc1 <- roc(df_for_plot$x_val, df_for_plot$y_val)
plot(roc1,
    thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
    print.thres="best",
    print.auc=TRUE,
    main=drug)

