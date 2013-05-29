package dream;

import java.io.FileNotFoundException;
import java.io.IOException;

public class Dream {

    String DATA_PATH = "/home/Dream/";

    public void makeFilteredGeneListForICCorrelation(int min_cnt) throws FileNotFoundException, IOException { // 120927
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
//        String sig_gene_list = DATA_PATH + "q2/count_drugs_per_sigGene_filtered.txt";
//        String sig_gene_list = DATA_PATH + "q2/count_drugs_per_sigGene_sorted.txt";   // probe based
        String sig_gene_list = DATA_PATH + "q2/count_drugs_per_Gene_sorted.txt"; // gene based
//        String input_dir = DATA_PATH + "q2/drug/";    // probe_based
        String input_dir = DATA_PATH + "q2/drug/gene_based/";   // gene_based
//        String output_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/";  // probe_based
        String output_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/"; // gene_based

        dh.makeFilteredGeneListForICCorrelation(drug_name_list, sig_gene_list, input_dir, output_dir, min_cnt);
    }

    public void makeFilteredGeneListForICCorrelation_with_target_gene(int min_cnt) throws FileNotFoundException, IOException { // 120927
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
        String target_gene_list = DATA_PATH + "q2/drug_target_gene.txt";
        String sig_gene_list = DATA_PATH + "q2/count_drugs_per_Gene_sorted.txt"; // gene based
        String input_dir = DATA_PATH + "q2/drug/gene_based/";   // gene_based
        String output_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/"; // gene_based

        dh.makeFilteredGeneListForICCorrelation_with_target_gene(drug_name_list, sig_gene_list, target_gene_list, input_dir, output_dir, min_cnt);
    }

    public void makeICCorrelationMatrix() throws FileNotFoundException, IOException { // 120927
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
        String input_dir = DATA_PATH + "q2/drug/filtered/";
        String output = DATA_PATH + "q2/drug/filtered/correlation_matrix.txt";

        dh.makeICCorrelationMatrix(drug_name_list, input_dir, output);
    }

    public void findCorrelatedDrugPairs() throws FileNotFoundException, IOException { // 120927
        DataHandler dh = new DataHandler();
        String score_matrix = DATA_PATH + "q2/drug/filtered/correlation_matrix.txt";
        String output = DATA_PATH + "q2/drug/filtered/correlated_drug_pairs.txt";
        double thres = 0.5;

        dh.findCorrelatedDrugPairs(score_matrix, output, thres);
    }

    public void makeFiltered24SigGeneList(int min_cnt) throws FileNotFoundException, IOException { // 120927
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
        String sensitive_drug_list = DATA_PATH + "q2/sensitive_drug_name.txt";
        String sensitive_drug_dir = DATA_PATH + "q2/sensitive/";
        String input_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/";
        String output_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/sig/";

        dh.makeFiltered24SigGeneList(drug_name_list, sensitive_drug_list, sensitive_drug_dir, input_dir, output_dir);
    }

    public void makeFiltered24SigGeneList_with_target_gene(int min_cnt) throws FileNotFoundException, IOException { // 120929
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
        String target_gene_list = DATA_PATH + "q2/drug_target_gene.txt";
        String input_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/";
        String output_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/";

        dh.makeFiltered24SigGeneList_with_target_gene(drug_name_list, target_gene_list, input_dir, output_dir);
    }

    public void makeDrugComparisonStatus(int min_cnt) throws FileNotFoundException, IOException { // 120927
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
//        String input_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/sig/"; // probe based
//        String output_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/sig/combination/";
        String input_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/"; // gene based
        String output_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/combination/";

        dh.makeDrugComparisonStatus(drug_name_list, input_dir, output_dir);
    }

    public void makeGeneBasedData() throws FileNotFoundException, IOException { // 120927
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
        String probe_to_gene = DATA_PATH + "q2/drug/Vincristine.txt";
        String input_dir = DATA_PATH + "q2/drug/";
        String output_dir = DATA_PATH + "q2/drug/gene_based/";

        dh.makeGeneBasedData(drug_name_list, probe_to_gene, input_dir, output_dir);
    }

    public void makeDrugCombinationScore(int min_cnt) throws FileNotFoundException, IOException { // 120929
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
        String input_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/";
        String output = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/combination/drug_combination_score.txt";

        dh.makeDrugCombinationScore(drug_name_list, input_dir, output);
    }

    public void checkTransporterGenes(int min_cnt) throws FileNotFoundException, IOException { // 120929
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
        String tr_gene_list = DATA_PATH + "q2/tr_gene_list.txt";
        String input_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/";
        String output = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/tr_gene_combination.txt";

        dh.checkTransporterGenes(drug_name_list, tr_gene_list, input_dir, output);
    }

    public void makeFilteredGeneList_for_tr_gene(int min_cnt) throws FileNotFoundException, IOException { // 121001
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
        String tr_gene_list = DATA_PATH + "q2/drug_resistance.txt";
        String input_dir = DATA_PATH + "q2/drug/gene_based/";   // gene_based
        String output_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/"; // gene_based

        dh.makeFilteredGeneList_for_tr_gene(drug_name_list, tr_gene_list, input_dir, output_dir);
    }

    public void makeDrugComparisonStatus_with_tr_gene(int min_cnt) throws FileNotFoundException, IOException { // 121001
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
        String tr_gene_list = DATA_PATH + "q2/drug_resistance.txt";
        String input_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/"; // gene based
        String output_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/combination/";
//        String output_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/same/";

        dh.makeDrugComparisonStatus_with_tr_gene(drug_name_list, tr_gene_list, input_dir, output_dir);
    }

    public void makeDrugCombinationScore_with_tr_gene(int min_cnt) throws FileNotFoundException, IOException { // 120929
        DataHandler dh = new DataHandler();
        String drug_name_list = DATA_PATH + "q2/drug_name.txt";
        String input_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/";
        String tr_gene_input_dir = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/";
        String output = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/combination/drug_combination_score_with_tr_gene.txt";
        String output_sep = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/combination/drug_combination_score_with_tr_gene_antagonist_separate.txt";
        String output_same = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/same/drug_combination_score_with_same_gene.txt";
        String output_sqrt = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/combination/drug_combination_score_with_tr_gene_sqrt.txt";
        String output_reverse = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/combination/drug_combination_score_with_tr_gene_reverse.txt";
        String output_full = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/combination/drug_combination_score_with_tr_gene_full.txt";
        String output_t_diff = DATA_PATH + "q2/drug/filtered/min_cnt_" + min_cnt + "/gene_based/sig/combination/drug_combination_score_with_tr_gene_t_diff.txt";

        dh.makeDrugCombinationScore_with_tr_gene(drug_name_list, input_dir, tr_gene_input_dir, output);
//        dh.makeDrugCombinationScore_with_tr_gene_antagonist_separate(drug_name_list, input_dir, tr_gene_input_dir, output_sep);
//        dh.makeDrugCombinationScore_with_tr_gene_between_same_gene(drug_name_list, input_dir, tr_gene_input_dir, output_same);
//        dh.makeDrugCombinationScore_with_tr_gene(drug_name_list, input_dir, tr_gene_input_dir, output_reverse);
//        dh.makeDrugCombinationScore_with_tr_gene(drug_name_list, input_dir, tr_gene_input_dir, output_full);
//        dh.makeDrugCombinationScore_with_tr_gene(drug_name_list, input_dir, tr_gene_input_dir, output_t_diff);
    }

    public static void main(String[] args) throws FileNotFoundException, IOException {

        Dream main = new Dream();
       // main.makeFilteredGeneListForICCorrelation(1);
       // main.makeICCorrelationMatrix();
       // main.findCorrelatedDrugPairs();

       // main.makeFilteredGeneListForICCorrelation(2);
       // main.makeFiltered24SigGeneList(2);

       // main.makeFilteredGeneListForICCorrelation(3);
       // main.makeFiltered24SigGeneList(3);

       // main.makeGeneBasedData();
       // main.makeDrugComparisonStatus(3);


       main.makeFilteredGeneListForICCorrelation_with_target_gene(3);
       main.makeFiltered24SigGeneList_with_target_gene(3);
       // main.makeDrugComparisonStatus(3);

       // main.makeDrugCombinationScore(3);

       // main.checkTransporterGenes(3);

       main.makeFilteredGeneList_for_tr_gene(3);
        main.makeDrugComparisonStatus_with_tr_gene(3);
        main.makeDrugCombinationScore_with_tr_gene(3);


    }
}
