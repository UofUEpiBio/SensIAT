####################################################
######  Function to generate times of a given post-baseline assessment
######  Uses Ogata's Thinning Algorithm to generate times that will have
######    a prescribed intensity -- in this case, the estimated intensity function
######    from the ARC data

######  Inputs are the Visit Number, vectors of subjects' previous assessment times
######    & outcomes, a value lambda_star for the thinning, and
######    baseline intensity vector lambda0_t and fitted cox model for desired
######    intensity function
######  Outputs a data frame with the assessment time, event indicator,
######     previous time & outcome, lag time, and visit number for each subject

######  Lambda_star should be a value > the maximum for the prescribed intensity
Times_gen_fn=function(Start_vec=rep(0,N),
                      Visit=1,
                      Outcomes_vec,
                      lambda0_t,
                      lambda_star,
                      cox_model){

    Prev_outcome_df=data.frame(Prev_outcome=c(Outcomes_vec),
                               Visit_number = rep(Visit, length(Outcomes_vec)))

    predict_vec=predict(cox_model, newdata=Prev_outcome_df, type="lp",reference="zero")
    # "lp" means linear predictor
    # A reference of "zero" causes no centering to be done

    Times_vec=rep(NA,N)
    Event_vec=rep(1,N)

    for(i in 1:N){

        start=Start_vec[i]

        done=0
        while(done==0){

            ####  no kth visit if they had no (k-1)st visit
            if(Start_vec[i]==End+Visit-5){

                Times_vec[i]=End+Visit-4
                done=1

            }else{

                #####  generate a potential assessment time Tstar
                S=rexp(1,rate=lambda_star)
                Tstar=ceiling(start+S)

                #####  if the candidate assessment time is too late, end with no additional
                #####     assessment time for that subject
                if(Tstar >=  End+Visit-5){

                    Times_vec[i]=End+Visit-4
                    done=1

                }else{

                    ######  if the candidate assessment time is not too late, accept it with
                    ######    a probability based on the the value of the intent
                    lambda=lambda0_t[Tstar]*exp(predict_vec[i])

                    U=runif(1)
                    accept=1*(U < lambda/lambda_star)

                    if(accept==1){

                        Times_vec[i]=Tstar
                        done=1

                    }else{ start=Tstar }

                }
            }
        }

        #####  Mark as having no new assessment time if applicable
        if( Times_vec[i]==(End+Visit-4) ){ Event_vec[i]=0 }

    }

    df=data.frame(id=1:N,
                  time=Times_vec,
                  event=Event_vec,
                  Prev_outcome=Outcomes_vec,
                  prev_time=Start_vec,
                  Lag_time=Times_vec-Start_vec,
                  visit=Visit)

    return(df)

}
###############################
