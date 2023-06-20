######
#Main Function
######

#Run this to generate all the tables from the FAST Paper.
#Set working directory to the root of the FAST Sims Folder
#The results from each simulation will be stored as .csv files in the results folder.

source("fast_script.R")
source("comparative_methods.R")
source("./Power Comparative/Order1_Comparison.R")
source("./Power Comparative/Order2_Comparison.R")
source("./Power Comparative/Order3_Comparison.R")
source("./Power Comparative/Order4_Comparison.R")
source("./Power Comparative/All_Scenarios_Partial_Anomaly_Period_Sequential.R")

source("./Power & Delay/Anomaly_At_Start.R")
source("./Power & Delay/Anomaly_At_End.R")
source("./Power & Delay/All_Scenarios_Order1.R")
source("./Power & Delay/All_Scenarios_Order2.R")
source("./Power & Delay/All_Scenarios_Order3.R")
source("./Power & Delay/All_Scenarios_Order4.R")
source("./Power & Delay/Contaminated_Training_Data.R")

source("./Choice of ODE Order/Varying_m_Order1.R")
source("./Choice of ODE Order/Varying_m_Order2.R")
source("./Choice of ODE Order/Varying_m_Order3.R")
source("./Choice of ODE Order/Varying_m_Order4.R")
source("./Choice of ODE Order/BIC_Consistency.R")
#store results as csv files

write.csv(result_df_order1_comparison, "./Results/Table12.csv")
write.csv(result_df_order2_comparison, "./Results/Table13.csv")
write.csv(result_df_order3_comparison, "./Results/Table14.csv")
write.csv(result_df_order4_comparison, "./Results/Table15.csv")
write.csv(result_df_partial_anomaly_sequential, "./Results/Review_Reply_Table1.csv")
write.csv(result_df_order2_comparison, "./Results/Review_Reply_Table2.csv")

write.csv(result_df_anomaly_at_start, "./Results/Table2_At_Start.csv")
write.csv(result_df_anomaly_at_end, "./Results/Table2_At_End.csv")
write.csv(result_df_order1, "./Results/Table6.csv")
write.csv(result_df_order2, "./Results/Table1.csv")
write.csv(result_df_order3, "./Results/Table7.csv")
write.csv(result_df_order4, "./Results/Table8.csv")
write.csv(result_df_contaminated, "./Results/Table3.csv")

write.csv(result_df_varying_m_order1, "./Results/Table9.csv")
write.csv(result_df_varying_m_order2, "./Results/Table4.csv")
write.csv(result_df_varying_m_order3, "./Results/Table10.csv")
write.csv(result_df_varying_m_order4, "./Results/Table11.csv")
write.csv(result_df_bic, "./Results/Table5.csv")

