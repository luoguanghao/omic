



function(nn, sample_ls) {

    ### divide
    # nn = 0.8 ### change!!
    sub<-sample(1:length(sample_ls),round(length(sample_ls)*nn))
    t_sample = as.vector(sample_ls[sub])
    v_sample = as.vector(sample_ls[-sub])

    return(list( t_sample_ls=t_sample, v_sample_ls=v_sample))


}
























