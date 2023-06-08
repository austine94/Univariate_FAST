######
#Main Function
######

#Run this to generate all the tables from the FAST Paper. These will be stored in the results folder.

source("fast_script.R")
source("comparative_methods.R")
source("./Power Comparative/All_Scenarios_Partial_Anomaly_Period.R")
source("./Power Comparative/All_Scenarios_Partial_Anomaly_Period_Sequential.R")
source("./Power Comparative/All_Scenarios_Full_Anomaly_Period.R")
source("./Power Comparative/All_Scenarios_Full_Anomaly_Period_Sequential.R")

source("./Power & Delay/Anomaly_At_End.R")
source("./Power & Delay/Varying_Anomaly_Length.R")
source("./Power & Delay/All_Scenarios_Order1.R")
source("./Power & Delay/All_Scenarios_Order2.R")
source("./Power & Delay/All_Scenarios_Order3.R")
source("./Power & Delay/All_Scenarios_Order4.R")


source("./Choice of ODE Order/Varying_m.R")
source("./Choice of ODE Order/BIC_Consistency.R")

#store results as csv files

write.csv(result_df_full_anomaly, "./Results/full_period_anomaly_results.csv")
write.csv(result_df_partial_anomaly, "./Results/partial_period_anomaly_results.csv")
write.csv(result_df_full_anomaly_sequential, "./Results/full_period_anomaly_sequential_results.csv")
write.csv(result_df_partial_anomaly_sequential, "./Results/partial_period_anomaly_sequentual_results.csv")

write.csv(result_df_anomaly_at_end, "./Results/anomaly_at_end_results.csv")
write.csv(varying_s_df, "./Results/varying_anomaly_length.csv")
write.csv(result_df_order1, "./Results/order1.csv")
write.csv(result_df_order2, "./Results/order2.csv")
write.csv(result_df_order3, "./Results/order3.csv")
write.csv(result_df_order4, "./Results/order4.csv")

write.csv(result_df_bic, "./Results/BIC_Consistency.csv")
write.csv(result_df_varying_m, "./Results/Varying_m.csv")
