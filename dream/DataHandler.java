package dream;

import java.io.*;
import java.util.*;

public class DataHandler {

    public void makeFilteredGeneListForICCorrelation(String drug_name_list, String sig_gene_list, String input_dir, String output_dir, int min_cnt) throws FileNotFoundException, IOException {   // 120927
        FileReader fr = new FileReader(sig_gene_list);
        BufferedReader br = new BufferedReader(fr);
        Hashtable sig_gene_table = new Hashtable(10000);

        String temp = br.readLine();
        String[] s;
        while (temp != null) {
            s = temp.split("\t");
            if (Integer.parseInt(s[1]) >= min_cnt) {
                sig_gene_table.put(s[0], "");
            }
            temp = br.readLine();
        }
        br.close();
        fr.close();

        fr = new FileReader(drug_name_list);
        br = new BufferedReader(fr);
        temp = br.readLine();
        while (temp != null) {
//            for (int i = 0; i < 3; i++) { // probe_based
//                int time = 6;
//                if (i == 1) {
//                    time = 12;
//                } else if (i == 2) {
//                    time = 24;
//                }
            int time = 24;  // for gene_based
//                String raw_input = input_dir + temp + "_" + time + "_IC20_ttest_fdr_combine.txt"; // probe_based
//                String filtered_output = output_dir + temp + "_" + time + "_IC20_ttest_fdr_combine_filtered.txt";
            String raw_input = input_dir + temp + "_" + time + "_IC20_gene_based_ttest_fdr.txt";   // gene_based
            String filtered_output = output_dir + temp + "_" + time + "_IC20_gene_based_ttest_fdr_combine_filtered.txt";

            FileReader ffr = new FileReader(raw_input);
            BufferedReader bbr = new BufferedReader(ffr);

            FileWriter fw = new FileWriter(filtered_output);
            BufferedWriter bw = new BufferedWriter(fw);

            String t = bbr.readLine();  // header
            t = t.replaceAll("\"", "");
            bw.write(t + "\n");   // header
            t = bbr.readLine();
            while (t != null) {
                t = t.replaceAll("\"", "");
                s = t.split("\t");
                if (sig_gene_table.containsKey(s[0])) {
                    bw.write(t + "\n");
                }
                t = bbr.readLine();
            }

            bbr.close();
            ffr.close();

            bw.flush();
            bw.close();
            fw.close();
//            }
            temp = br.readLine();
        }
        br.close();
        fr.close();
    }

    public void makeFilteredGeneListForICCorrelation_with_target_gene(String drug_name_list, String sig_gene_list, String target_gene_list, String input_dir, String output_dir, int min_cnt) throws FileNotFoundException, IOException {   // 120927
        FileReader fr = new FileReader(sig_gene_list);
        BufferedReader br = new BufferedReader(fr);
        Hashtable sig_gene_table = new Hashtable(10000);

        String temp = br.readLine();
        String[] s;
        while (temp != null) {
            s = temp.split("\t");
            if (Integer.parseInt(s[1]) >= min_cnt) {
                sig_gene_table.put(s[0], "");
            }
            temp = br.readLine();
        }
        br.close();
        fr.close();

        fr = new FileReader(target_gene_list);
        br = new BufferedReader(fr);

        temp = br.readLine(); // header
        temp = br.readLine();
        while (temp != null) {
            temp = temp.replaceAll("\"", "");
            s = temp.split("\t");
            if (s.length > 1) {
                String[] g = s[1].split(",");
                if (g.length > 0) {
                    for (int i = 0; i < g.length; i++) {
                        g[i] = g[i].trim();
                        sig_gene_table.put(g[i], "");
                    }
                }
            }
            temp = br.readLine();
        }
        br.close();
        fr.close();

        fr = new FileReader(drug_name_list);
        br = new BufferedReader(fr);
        temp = br.readLine();
        while (temp != null) {
//            for (int i = 0; i < 3; i++) { // probe_based
//                int time = 6;
//                if (i == 1) {
//                    time = 12;
//                } else if (i == 2) {
//                    time = 24;
//                }
            int time = 24;  // for gene_based
//                String raw_input = input_dir + temp + "_" + time + "_IC20_ttest_fdr_combine.txt"; // probe_based
//                String filtered_output = output_dir + temp + "_" + time + "_IC20_ttest_fdr_combine_filtered.txt";
            String raw_input = input_dir + temp + "_" + time + "_IC20_gene_based_ttest_fdr.txt";   // gene_based
            String filtered_output = output_dir + temp + "_" + time + "_IC20_gene_based_ttest_fdr_combine_filtered.txt";

            FileReader ffr = new FileReader(raw_input);
            BufferedReader bbr = new BufferedReader(ffr);

            FileWriter fw = new FileWriter(filtered_output);
            BufferedWriter bw = new BufferedWriter(fw);

            String t = bbr.readLine();  // header
            t = t.replaceAll("\"", "");
            bw.write(t + "\n");   // header
            t = bbr.readLine();
            while (t != null) {
                t = t.replaceAll("\"", "");
                s = t.split("\t");
                if (sig_gene_table.containsKey(s[0])) {
                    bw.write(t + "\n");
                }
                t = bbr.readLine();
            }

            bbr.close();
            ffr.close();

            bw.flush();
            bw.close();
            fw.close();
//            }
            temp = br.readLine();
        }
        br.close();
        fr.close();
    }

    public void makeICCorrelationMatrix(String drug_name_list, String input_dir, String output) throws FileNotFoundException, IOException {   // 120927
        FileReader fr = new FileReader(drug_name_list);
        BufferedReader br = new BufferedReader(fr);

        LinkedList drug_list = new LinkedList();
        LinkedList name_list = new LinkedList();

        int drug_cnt = 0;

        String temp = br.readLine();
        String[] s;
        System.out.println("Making t_statics vectors ...");
        while (temp != null) {
            drug_cnt++;
            LinkedList time_list = new LinkedList();
            for (int i = 0; i < 3; i++) {
                int time = 6;
                if (i == 1) {
                    time = 12;
                } else if (i == 2) {
                    time = 24;
                }
                String raw_input = input_dir + temp + "_" + time + "_IC20_ttest_fdr_combine_filtered.txt";

                FileReader ffr = new FileReader(raw_input);
                BufferedReader bbr = new BufferedReader(ffr);

                LinkedList t_statistics = new LinkedList();

                String t = bbr.readLine();  // header
                t = bbr.readLine();
                while (t != null) {
                    t = t.replaceAll("\"", "");
                    s = t.split("\t");
                    t_statistics.add(s[1]);
                    t = bbr.readLine();
                }
                bbr.close();
                ffr.close();

                time_list.add(t_statistics);
            }

            drug_list.add(time_list);
            name_list.add(temp);

            temp = br.readLine();
        }
        br.close();
        fr.close();

        FileWriter fw = new FileWriter(output);
        BufferedWriter bw = new BufferedWriter(fw);
        String header = "Drug_list";
        for (int i = 0; i < name_list.size(); i++) {
            String drug = (String) name_list.get(i);
            int num = 6;
            for (int j = 0; j < 3; j++) {
                header += "\t" + drug + "_" + num;
                num = num * 2;
            }
        }
        bw.write(header + "\n");

        double[][] score_matrix = new double[drug_cnt * 3][drug_cnt * 3];
        System.out.println("Making Scoring matrix ...");

        for (int i = 0; i < drug_cnt; i++) {
            LinkedList drug_i = (LinkedList) drug_list.get(i);
            int current_row = i * 3;
            System.out.print("- drug_" + (i + 1));
            for (int j = 0; j < drug_cnt; j++) {
                System.out.print(" - drug_" + (j + 1));
                LinkedList drug_j = (LinkedList) drug_list.get(j);
                int current_col = j * 3;
                if (i == j) {
                    for (int d_i = 0; d_i < 3; d_i++) {
                        for (int d_j = 0; d_j < 3; d_j++) {
                            score_matrix[current_row + d_i][current_col + d_j] = 0;
                        }
                    }
                } else {
                    for (int d_i = 0; d_i < 3; d_i++) {
                        LinkedList t_i = (LinkedList) drug_i.get(d_i);
                        for (int d_j = 0; d_j < 3; d_j++) {
                            LinkedList t_j = (LinkedList) drug_j.get(d_j);
                            double cor = getCorrelation(t_i, t_j);
                            score_matrix[current_row + d_i][current_col + d_j] = cor;
                        }
                    }
                }
            }
            System.out.println();
        }

        System.out.println("Making output file ...");

        for (int i = 0; i < drug_cnt; i++) {
            String drug = (String) name_list.get(i);
            int num = 6;
            int current_row = i * 3;
            for (int j = 0; j < 3; j++) {
                String result = drug + "_" + num;
                num = num * 2;

                for (int k = 0; k < drug_cnt * 3; k++) {
                    result += "\t" + score_matrix[current_row + j][k];
                }
                bw.write(result + "\n");
            }
        }

        bw.flush();
        bw.close();
        fw.close();
    }

    public void findCorrelatedDrugPairs(String score_matrix, String output, double thres) throws FileNotFoundException, IOException {   // 120927
        FileReader fr = new FileReader(score_matrix);
        BufferedReader br = new BufferedReader(fr);

        FileWriter fw = new FileWriter(output);
        BufferedWriter bw = new BufferedWriter(fw);

        LinkedList drug_list = new LinkedList();

        String temp = br.readLine();    // header
        String[] s;
        s = temp.split("\t");
        for (int i = 0; i < s.length; i++) {
            drug_list.add(s[i]);
        }
        temp = br.readLine();
        int current_row = 1;
        while (temp != null) {
            s = temp.split("\t");
            for (int i = current_row + 1; i < s.length; i++) {
                if (Math.abs(Double.parseDouble(s[i])) >= thres) {
                    String src = (String) drug_list.get(current_row);
                    String tgt = (String) drug_list.get(i);
                    bw.write(src + "\t" + tgt + "\t" + s[i] + "\n");
                }
            }
            temp = br.readLine();
            current_row++;
        }
        br.close();
        fr.close();

        bw.flush();
        bw.close();
        fw.close();
    }

    public void test() {
        double[] x = {2, 3, 4, 5, 6, 8, 10, 11.53542};
        double[] y = {21.05, 23.51, 24.23, 27.71, 30.86, 45.85, 52.12, 55.98};

        LinkedList x_list = new LinkedList();
        LinkedList y_list = new LinkedList();
        for (int i = 0; i < x.length; i++) {
            x_list.add(Double.toString(x[i]));
            y_list.add(Double.toString(y[i]));
        }

        System.out.println(getCorrelation(x_list, y_list));
    }

    public double getMean(LinkedList x) {
        double sum = 0;
        for (int i = 0; i < x.size(); i++) {
            sum += Double.parseDouble((String) x.get(i));
        }
        return sum / x.size();
    }

    public double getStandardDeviation(LinkedList x) {
        double mean = getMean(x);
        double sumOfSquaredDeviations = 0;

        for (int i = 0; i < x.size(); i++) {
            sumOfSquaredDeviations += Math.pow(Double.parseDouble((String) x.get(i)) - mean, 2);
        }
        return Math.sqrt(sumOfSquaredDeviations / (x.size() - 1));
    }

    public double getCovariance(LinkedList x, LinkedList y) {
        double result = 0;
        double xmean = getMean(x);
        double ymean = getMean(y);

        for (int i = 0; i < x.size(); i++) {
            result += (Double.parseDouble((String) x.get(i)) - xmean) * (Double.parseDouble((String) y.get(i)) - ymean);
        }
        result /= x.size() - 1;

        return result;
    }

    public double getCorrelation(LinkedList x, LinkedList y) {
        double xStdDev = getStandardDeviation(x);
        double yStdDev = getStandardDeviation(y);

        return getCovariance(x, y) / (xStdDev * yStdDev);
    }

    public void makeFiltered24SigGeneList(String drug_name_list, String sensitive_drug_list, String sensitive_drug_dir, String input_dir, String output_dir) throws FileNotFoundException, IOException {   // 120927
        FileReader fr = new FileReader(sensitive_drug_list);
        BufferedReader br = new BufferedReader(fr);
        Hashtable sensitive_drug_table = new Hashtable(5);

        String temp = br.readLine();
        while (temp != null) {
            sensitive_drug_table.put(temp, "");
            temp = br.readLine();
        }
        br.close();
        fr.close();

        fr = new FileReader(drug_name_list);
        br = new BufferedReader(fr);
        temp = br.readLine();
        String[] s;
        while (temp != null) {
            int time = 24;
            String raw_input = input_dir + temp + "_" + time + "_IC20_ttest_fdr_combine_filtered.txt";
            String sig_filtered_output = output_dir + temp + "_" + time + "_IC20_sig_filtered.txt";
            Hashtable sensitive_gene_table = new Hashtable(20000);

            if (sensitive_drug_table.containsKey(temp)) {
                FileReader ffr = new FileReader(sensitive_drug_dir + temp + "_sensitive_gene.txt");
                BufferedReader bbr = new BufferedReader(ffr);
                String t = bbr.readLine();
                while (t != null) {
                    t = t.replaceAll("\"", "");
                    sensitive_gene_table.put(t, "");
                    t = bbr.readLine();
                }
                bbr.close();
                ffr.close();
            }

            FileReader ffr = new FileReader(raw_input);
            BufferedReader bbr = new BufferedReader(ffr);

            FileWriter fw = new FileWriter(sig_filtered_output);
            BufferedWriter bw = new BufferedWriter(fw);

            String t = bbr.readLine();  // header
            String header = "sig_gene\tt-score\tdirection\tisSensitive";
            bw.write(header + "\n");   // header
            t = bbr.readLine();
            while (t != null) {
                t = t.replaceAll("\"", "");
                s = t.split("\t");
                if (Double.parseDouble(s[4]) < 0.1) {
                    boolean isSensitive = false;
                    if (sensitive_gene_table.containsKey(s[0])) {
                        isSensitive = true;
                    }
                    String direction = "+";
                    if (Double.parseDouble(s[1]) < 0) {
                        direction = "-";
                    }
                    String result = s[0] + "\t" + s[1] + "\t" + direction + "\t" + isSensitive;
                    bw.write(result + "\n");
                }
                t = bbr.readLine();
            }

            bbr.close();
            ffr.close();

            bw.flush();
            bw.close();
            fw.close();

            sensitive_gene_table.clear();

            temp = br.readLine();
        }
        br.close();
        fr.close();
    }

    public void makeFiltered24SigGeneList_with_target_gene(String drug_name_list, String target_gene_list, String input_dir, String output_dir) throws FileNotFoundException, IOException {   // 120927
        FileReader fr = new FileReader(target_gene_list);
        BufferedReader br = new BufferedReader(fr);

        Hashtable target_table = new Hashtable(20);

        String temp = br.readLine();    // header
        temp = br.readLine();
        String[] s;
        while (temp != null) {
            temp = temp.replaceAll("\"", "");
            s = temp.split("\t");
            Hashtable target_gene_table = new Hashtable(10);
            if (s.length > 1) {
                String[] g = s[1].split(",");
                if (g.length > 0) {
                    for (int i = 0; i < g.length; i++) {
                        g[i] = g[i].trim();
                        target_gene_table.put(g[i], "");
                    }
                }
            }
            target_table.put(s[0], target_gene_table);

            temp = br.readLine();
        }
        br.close();
        fr.close();

        fr = new FileReader(drug_name_list);
        br = new BufferedReader(fr);
        temp = br.readLine();
        while (temp != null) {
            int time = 24;
            String raw_input = input_dir + temp + "_" + time + "_IC20_gene_based_ttest_fdr_combine_filtered.txt";
            String sig_filtered_output = output_dir + temp + "_" + time + "_IC20_gene_based_sig_filtered.txt";

            FileReader ffr = new FileReader(raw_input);
            BufferedReader bbr = new BufferedReader(ffr);

            FileWriter fw = new FileWriter(sig_filtered_output);
            BufferedWriter bw = new BufferedWriter(fw);

            String t = bbr.readLine();  // header
            String header = "sig_gene\tt-score\tdirection\tisSensitive";
            bw.write(header + "\n");   // header
            t = bbr.readLine();
            Hashtable target_gene_table = (Hashtable) target_table.get(temp);
            while (t != null) {
                t = t.replaceAll("\"", "");
                s = t.split("\t");

                if (target_gene_table.containsKey(s[0])) {
                    boolean isSensitive = false;
                    String direction = "+";
                    if (Double.parseDouble(s[1]) < 0) {
                        direction = "-";
                    }
                    String result = s[0] + "\t" + s[1] + "\t" + direction + "\t" + isSensitive;
                    bw.write(result + "\n");
                } else if (Double.parseDouble(s[4]) < 0.1) {
                    boolean isSensitive = false;
                    String direction = "+";
                    if (Double.parseDouble(s[1]) < 0) {
                        direction = "-";
                    }
                    String result = s[0] + "\t" + s[1] + "\t" + direction + "\t" + isSensitive;
                    bw.write(result + "\n");
                }
                t = bbr.readLine();
            }

            bbr.close();
            ffr.close();

            bw.flush();
            bw.close();
            fw.close();

            temp = br.readLine();
        }
        br.close();
        fr.close();
    }

    public void makeDrugComparisonStatus(String drug_name_list, String input_dir, String output_dir) throws FileNotFoundException, IOException {   // 120927
        FileReader fr = new FileReader(drug_name_list);
        BufferedReader br = new BufferedReader(fr);

        LinkedList drug_list = new LinkedList();

        String temp = br.readLine();
        while (temp != null) {
            drug_list.add(temp);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        for (int i = 0; i < drug_list.size() - 1; i++) {
            for (int j = i + 1; j < drug_list.size(); j++) {
                makeComparisonOfDrugPair((String) drug_list.get(i), (String) drug_list.get(j), input_dir, output_dir);
            }
        }

    }

    public void makeComparisonOfDrugPair(String drug_x, String drug_y, String input_dir, String output_dir) throws FileNotFoundException, IOException {   // 120927
//        FileReader fr = new FileReader(input_dir + drug_x + "_24_IC20_sig_filtered.txt");   // probe_based
        FileReader fr = new FileReader(input_dir + drug_x + "_24_IC20_gene_based_sig_filtered.txt");   // gene_based
        BufferedReader br = new BufferedReader(fr);

        LinkedList only_x = new LinkedList();
        LinkedList only_y = new LinkedList();
        LinkedList common_same_d_sensitive = new LinkedList();
        LinkedList common_same_d_non_sensitive = new LinkedList();
        LinkedList common_opposite_d_sensitive = new LinkedList();
        LinkedList common_opposite_d_non_sensitive = new LinkedList();

        FileWriter fw = new FileWriter(output_dir + drug_x + "_" + drug_y + "_combination_stats.txt");
        BufferedWriter bw = new BufferedWriter(fw);

        Hashtable x_table = new Hashtable(8000);

        String temp = br.readLine();    // header
        temp = br.readLine();
        String[] s;
        while (temp != null) {
            s = temp.split("\t");
            x_table.put(s[0], temp);

            temp = br.readLine();
        }
        br.close();
        fr.close();

//        fr = new FileReader(input_dir + drug_y + "_24_IC20_sig_filtered.txt");  // probe_based
        fr = new FileReader(input_dir + drug_y + "_24_IC20_gene_based_sig_filtered.txt");  // probe_based
        br = new BufferedReader(fr);
        temp = br.readLine();   // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            String data = (String) x_table.remove(s[0]);
            if (data == null) {
                only_y.add(s[0]);
            } else {
                String[] s_x = data.split("\t");
                String d_x = s_x[2];
                String d_y = s[2];
                String isSensitive_x = s_x[3];
                String isSensitive_y = s[3];
                if (isSensitive_x.equals("true") || isSensitive_y.equals("true")) {
                    if (d_x.equals(d_y)) {
                        common_same_d_sensitive.add(s[0]);
                    } else {
                        common_opposite_d_sensitive.add(s[0]);
                    }
                } else {
                    if (d_x.equals(d_y)) {
                        common_same_d_non_sensitive.add(s[0]);
                    } else {
                        common_opposite_d_non_sensitive.add(s[0]);
                    }
                }
            }

            temp = br.readLine();
        }
        br.close();
        fr.close();

        Set set = x_table.keySet();
        LinkedList keylist = new LinkedList(set);
        Iterator keys = keylist.iterator();
        while (keys.hasNext()) {
            String key = (String) keys.next();
            String x = (String) x_table.get(key);
            s = x.split("\t");
            only_x.add(s[0]);
        }

        String result = "";

        bw.write("> Genes only from " + drug_x + "\n");
        if (only_x.size() == 0) {
            result = "0";
        } else {
            result = only_x.size() + "\t" + only_x.get(0);
            for (int i = 1; i < only_x.size(); i++) {
                result += "," + only_x.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write("> Genes only from " + drug_y + "\n");
        if (only_y.size() == 0) {
            result = "0";
        } else {
            result = only_y.size() + "\t" + only_y.get(0);
            for (int i = 1; i < only_y.size(); i++) {
                result += "," + only_y.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write("> Common genes from both drugs\n");
        bw.write(">> same direction, sensitive genes\n");
        if (common_same_d_sensitive.size() == 0) {
            result = "0";
        } else {
            result = common_same_d_sensitive.size() + "\t" + common_same_d_sensitive.get(0);
            for (int i = 1; i < common_same_d_sensitive.size(); i++) {
                result += "," + common_same_d_sensitive.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write(">> same direction, non_sensitive genes\n");
        if (common_same_d_non_sensitive.size() == 0) {
            result = "0";
        } else {
            result = common_same_d_non_sensitive.size() + "\t" + common_same_d_non_sensitive.get(0);
            for (int i = 1; i < common_same_d_non_sensitive.size(); i++) {
                result += "," + common_same_d_non_sensitive.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write(">> opposite direction, sensitive genes\n");
        if (common_opposite_d_sensitive.size() == 0) {
            result = "0";
        } else {
            result = common_opposite_d_sensitive.size() + "\t" + common_opposite_d_sensitive.get(0);
            for (int i = 1; i < common_opposite_d_sensitive.size(); i++) {
                result += "," + common_opposite_d_sensitive.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write(">> opposite direction, non_sensitive genes\n");
        if (common_opposite_d_non_sensitive.size() == 0) {
            result = "0";
        } else {
            result = common_opposite_d_non_sensitive.size() + "\t" + common_opposite_d_non_sensitive.get(0);
            for (int i = 1; i < common_opposite_d_non_sensitive.size(); i++) {
                result += "," + common_opposite_d_non_sensitive.get(i);
            }
        }
        bw.write(result + "\n");

        bw.flush();
        bw.close();
        fw.close();
    }

    public void makeGeneBasedData(String drug_name_list, String probe_to_gene, String input_dir, String output_dir) throws FileNotFoundException, IOException {   // 120927
        FileReader fr = new FileReader(probe_to_gene);
        BufferedReader br = new BufferedReader(fr);
        Hashtable probe_to_gene_table = new Hashtable(50000);

        String temp = br.readLine();
        temp = br.readLine();
        temp = br.readLine();
        temp = br.readLine();    // header
        String[] s;
        while (temp != null) {
            s = temp.split("\t");
            probe_to_gene_table.put(s[0], s[1]);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        fr = new FileReader(drug_name_list);
        br = new BufferedReader(fr);
        temp = br.readLine();
        while (temp != null) {
            int time = 24;
            String raw_input = input_dir + temp + "_" + time + "_IC20.txt";
            String gene_output = output_dir + temp + "_" + time + "_IC20_gene_based.txt";

            FileReader ffr = new FileReader(raw_input);
            BufferedReader bbr = new BufferedReader(ffr);

            FileWriter fw = new FileWriter(gene_output);
            BufferedWriter bw = new BufferedWriter(fw);

            Hashtable data_table = new Hashtable(50000);

            String t = bbr.readLine();  // header
            t = t.replaceAll("\"", "");
            bw.write(t + "\n");   // header
            t = bbr.readLine();
            while (t != null) {
                t = t.replaceAll("\"", "");
                s = t.split("\t");

                String gene = (String) probe_to_gene_table.get(s[0]);
                boolean isSet = false;
                String[] gene_set = null;
                if (gene.indexOf("///") != -1) {
                    gene_set = gene.split(" /// ");
                    isSet = true;
                }

                if (!isSet) {
                    String data = (String) data_table.get(gene);
                    if (data == null) {
                        data = t.replace(s[0], gene);
                        data_table.put(gene, data);
                    } else {
                        String[] d = data.split("\t");
                        for (int i = 1; i < d.length; i++) {
                            if (Double.parseDouble(d[i]) < Double.parseDouble(s[i])) {
                                d[i] = s[i];
                            }
                        }
                        data = d[0];
                        for (int i = 1; i < d.length; i++) {
                            data += "\t" + d[i];
                        }
                        data_table.put(gene, data);
                    }
                } else {
                    for (int i = 0; i < gene_set.length; i++) {
                        String data = (String) data_table.get(gene_set[i]);
                        if (data == null) {
                            data = t.replace(s[0], gene_set[i]);
                            data_table.put(gene_set[i], data);
                        } else {
                            String[] d = data.split("\t");
                            for (int j = 1; j < d.length; j++) {
                                if (Double.parseDouble(d[j]) < Double.parseDouble(s[j])) {
                                    d[j] = s[j];
                                }
                            }
                            data = d[0];
                            for (int j = 1; j < d.length; j++) {
                                data += "\t" + d[j];
                            }
                            data_table.put(gene_set[i], data);
                        }
                    }
                }
                t = bbr.readLine();
            }

            bbr.close();
            ffr.close();

            Set set = data_table.keySet();
            LinkedList keylist = new LinkedList(set);
            Collections.sort(keylist);
            Iterator keys = keylist.iterator();
            while (keys.hasNext()) {
                String key = (String) keys.next();
                String result = (String) data_table.get(key);
                bw.write(result + "\n");
            }

            bw.flush();
            bw.close();
            fw.close();

            data_table.clear();

            temp = br.readLine();
        }
        br.close();
        fr.close();
    }

    public void makeDrugCombinationScore(String drug_name_list, String input_dir, String output) throws FileNotFoundException, IOException {   // 120927
        FileReader fr = new FileReader(drug_name_list);
        BufferedReader br = new BufferedReader(fr);

        FileWriter fw = new FileWriter(output);
        BufferedWriter bw = new BufferedWriter(fw);

        LinkedList drug_list = new LinkedList();

        String temp = br.readLine();
        while (temp != null) {
            drug_list.add(temp);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        for (int i = 0; i < drug_list.size() - 1; i++) {
            for (int j = i + 1; j < drug_list.size(); j++) {
                String combination = (String) drug_list.get(i) + "_" + (String) drug_list.get(j);
                String filename = input_dir + "combination/" + combination + "_combination_stats.txt";
                DrugCombination dc = new DrugCombination();
                dc.makeDrugCombinationProfile(filename);
                double score = dc.getCombinationScore(input_dir);
                bw.write(combination + "\t" + score + "\n");
            }
        }

        bw.flush();
        bw.close();
        fw.close();
    }

    public void checkTransporterGenes(String drug_name_list, String tr_gene_list, String input_dir, String output) throws FileNotFoundException, IOException {   // 120927
        FileReader fr = new FileReader(tr_gene_list);
        BufferedReader br = new BufferedReader(fr);

        Hashtable tr_gene_table = new Hashtable(100);

        String temp = br.readLine();
        while (temp != null) {
            LinkedList drug_list = new LinkedList();
            tr_gene_table.put(temp, drug_list);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        fr = new FileReader(drug_name_list);
        br = new BufferedReader(fr);
        temp = br.readLine();
        String[] s;
        while (temp != null) {
            int time = 24;
            String raw_input = input_dir + temp + "_" + time + "_IC20_gene_based_sig_filtered.txt";

            FileReader ffr = new FileReader(raw_input);
            BufferedReader bbr = new BufferedReader(ffr);

            String t = bbr.readLine();  // header
            t = bbr.readLine();
            while (t != null) {
                t = t.replaceAll("\"", "");
                s = t.split("\t");

                LinkedList drug_list = (LinkedList) tr_gene_table.get(s[0]);
                if (drug_list != null) {
                    if (s[2].equals("+")) {
                        drug_list.add(temp);
                    }
                }
                t = bbr.readLine();
            }
            bbr.close();
            ffr.close();

            temp = br.readLine();
        }
        br.close();
        fr.close();

        FileWriter fw = new FileWriter(output);
        BufferedWriter bw = new BufferedWriter(fw);

        Set set = tr_gene_table.keySet();
        LinkedList keylist = new LinkedList(set);
        Collections.sort(keylist);
        Iterator keys = keylist.iterator();
        while (keys.hasNext()) {
            String key = (String) keys.next();
            LinkedList drug_list = (LinkedList) tr_gene_table.get(key);
            String result = key;
            for (int i = 0; i < drug_list.size(); i++) {
                result += "\t" + drug_list.get(i);
            }
            bw.write(result + "\n");
        }

        bw.flush();
        bw.close();
        fw.close();
    }

    public void makeFilteredGeneList_for_tr_gene(String drug_name_list, String tr_gene_list, String input_dir, String output_dir) throws FileNotFoundException, IOException {   // 121001
        FileReader fr = new FileReader(tr_gene_list);
        BufferedReader br = new BufferedReader(fr);
        Hashtable tr_gene_table = new Hashtable(10000);

        String temp = br.readLine();    // header
        String[] s;
        temp = br.readLine();
        while (temp != null) {
            temp = temp.replaceAll("\"", "");
            s = temp.split("\t");
            if (s.length > 1) {
                String[] g = s[1].split(",");
                if (g.length > 0) {
                    for (int i = 0; i < g.length; i++) {
                        g[i] = g[i].trim();
                        tr_gene_table.put(g[i], "");
                    }
                }
                if (s.length > 2) {
                    g = s[2].split(",");
                    if (g.length > 0) {
                        for (int i = 0; i < g.length; i++) {
                            g[i] = g[i].trim();
                            tr_gene_table.put(g[i], "");
                        }
                    }
                }
            }
            temp = br.readLine();
        }
        br.close();
        fr.close();

        fr = new FileReader(drug_name_list);
        br = new BufferedReader(fr);
        temp = br.readLine();
        while (temp != null) {
//            for (int i = 0; i < 3; i++) { // probe_based
//                int time = 6;
//                if (i == 1) {
//                    time = 12;
//                } else if (i == 2) {
//                    time = 24;
//                }
            int time = 24;  // for gene_based
//                String raw_input = input_dir + temp + "_" + time + "_IC20_ttest_fdr_combine.txt"; // probe_based
//                String filtered_output = output_dir + temp + "_" + time + "_IC20_ttest_fdr_combine_filtered.txt";
            String raw_input = input_dir + temp + "_" + time + "_IC20_gene_based_ttest_fdr.txt";   // gene_based
            String filtered_output = output_dir + temp + "_" + time + "_IC20_gene_based_tr_gene_filtered.txt";

            FileReader ffr = new FileReader(raw_input);
            BufferedReader bbr = new BufferedReader(ffr);

            FileWriter fw = new FileWriter(filtered_output);
            BufferedWriter bw = new BufferedWriter(fw);

            String t = bbr.readLine();  // header
            t = t.replaceAll("\"", "");
            bw.write(t + "\n");   // header
            t = bbr.readLine();
            while (t != null) {
                t = t.replaceAll("\"", "");
                s = t.split("\t");
                if (tr_gene_table.containsKey(s[0])) {
                    bw.write(t + "\n");
                }
                t = bbr.readLine();
            }

            bbr.close();
            ffr.close();

            bw.flush();
            bw.close();
            fw.close();
//            }
            temp = br.readLine();
        }
        br.close();
        fr.close();
    }

    public void makeDrugComparisonStatus_with_tr_gene(String drug_name_list, String tr_gene_list, String input_dir, String output_dir) throws FileNotFoundException, IOException {   // 121001
        FileReader fr = new FileReader(drug_name_list);
        BufferedReader br = new BufferedReader(fr);

        LinkedList drug_list = new LinkedList();

        String temp = br.readLine();
        while (temp != null) {
            drug_list.add(temp);
            temp = br.readLine();
        }
        br.close();
        fr.close();        
        Collections.sort(drug_list);

        Hashtable inhibit_table = new Hashtable(20);
        Hashtable activate_table = new Hashtable(20);
        fr = new FileReader(tr_gene_list);
        br = new BufferedReader(fr);
        temp = br.readLine();   // header
        temp = br.readLine();
        String[] s;
        while (temp != null) {
            temp = temp.replaceAll("\"", "");
            s = temp.split("\t");
            LinkedList inhibit_list = new LinkedList();
            LinkedList activate_list = new LinkedList();
            if (s.length > 1) {
                String[] g = s[1].split(",");
                if (g.length > 0) {
                    for (int i = 0; i < g.length; i++) {
                        g[i] = g[i].trim();
                        inhibit_list.add(g[i]);
                    }
                }
                if (s.length > 2) {
                    g = s[2].split(",");
                    if (g.length > 0) {
                        for (int i = 0; i < g.length; i++) {
                            g[i] = g[i].trim();
                            activate_list.add(g[i]);
                        }
                    }
                }
            }
            s[0] = s[0].trim();
            inhibit_table.put(s[0], inhibit_list);
            activate_table.put(s[0], activate_list);
            temp = br.readLine();
        }
        br.close();
        fr.close();

//        for (int i = 0; i < drug_list.size(); i++) {
        for (int i = 0; i < drug_list.size() - 1; i++) {
            for (int j = i + 1; j < drug_list.size(); j++) {
                makeComparisonOfDrugPair_with_tr_gene((String) drug_list.get(i), (String) drug_list.get(j), inhibit_table, activate_table, input_dir, output_dir);
//                makeComparisonOfDrugPair_with_tr_gene((String) drug_list.get(i), (String) drug_list.get(i), inhibit_table, activate_table, input_dir, output_dir);
            }
        }

    }

    public void makeComparisonOfDrugPair_with_tr_gene(String drug_x, String drug_y, Hashtable inhibit_table, Hashtable activate_table, String input_dir, String output_dir) throws FileNotFoundException, IOException {   // 120927
//        FileReader fr = new FileReader(input_dir + drug_x + "_24_IC20_sig_filtered.txt");   // probe_based
        FileReader fr = new FileReader(input_dir + drug_x + "_24_IC20_gene_based_sig_filtered.txt");   // gene_based
        BufferedReader br = new BufferedReader(fr);

        LinkedList only_x = new LinkedList();
        LinkedList only_y = new LinkedList();
        LinkedList common_same_d_sensitive = new LinkedList();
        LinkedList common_same_d_non_sensitive = new LinkedList();
        LinkedList common_opposite_d_sensitive = new LinkedList();
        LinkedList common_opposite_d_non_sensitive = new LinkedList();
        LinkedList inhibit_x;
        LinkedList activate_x;
        LinkedList inhibit_y;
        LinkedList activate_y;

        FileWriter fw = new FileWriter(output_dir + drug_x + "_" + drug_y + "_combination_stats.txt");
        BufferedWriter bw = new BufferedWriter(fw);

        Hashtable x_table = new Hashtable(8000);

        String temp = br.readLine();    // header
        temp = br.readLine();
        String[] s;
        while (temp != null) {
            s = temp.split("\t");
            x_table.put(s[0], temp);

            temp = br.readLine();
        }
        br.close();
        fr.close();

//        fr = new FileReader(input_dir + drug_y + "_24_IC20_sig_filtered.txt");  // probe_based
        fr = new FileReader(input_dir + drug_y + "_24_IC20_gene_based_sig_filtered.txt");  // gene_based
        br = new BufferedReader(fr);
        temp = br.readLine();   // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            String data = (String) x_table.remove(s[0]);
            if (data == null) {
                only_y.add(s[0]);
            } else {
                String[] s_x = data.split("\t");
                String d_x = s_x[2];
                String d_y = s[2];
                String isSensitive_x = s_x[3];
                String isSensitive_y = s[3];
                if (isSensitive_x.equals("true") || isSensitive_y.equals("true")) {
                    if (d_x.equals(d_y)) {
                        common_same_d_sensitive.add(s[0]);
                    } else {
                        common_opposite_d_sensitive.add(s[0]);
                    }
                } else {
                    if (d_x.equals(d_y)) {
                        common_same_d_non_sensitive.add(s[0]);
                    } else {
                        common_opposite_d_non_sensitive.add(s[0]);
                    }
                }
            }

            temp = br.readLine();
        }
        br.close();
        fr.close();

        Set set = x_table.keySet();
        LinkedList keylist = new LinkedList(set);
        Iterator keys = keylist.iterator();
        while (keys.hasNext()) {
            String key = (String) keys.next();
            String x = (String) x_table.get(key);
            s = x.split("\t");
            only_x.add(s[0]);
        }

        String result = "";

        bw.write("> Genes only from " + drug_x + "\n");
        if (only_x.size() == 0) {
            result = "0";
        } else {
            result = only_x.size() + "\t" + only_x.get(0);
            for (int i = 1; i < only_x.size(); i++) {
                result += "," + only_x.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write("> Genes only from " + drug_y + "\n");
        if (only_y.size() == 0) {
            result = "0";
        } else {
            result = only_y.size() + "\t" + only_y.get(0);
            for (int i = 1; i < only_y.size(); i++) {
                result += "," + only_y.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write("> Inhibit transporter genes from " + drug_x + "\n");
        inhibit_x = (LinkedList) inhibit_table.get(drug_x);
        if (inhibit_x.size() == 0) {
            result = "0";
        } else {
            result = inhibit_x.size() + "\t" + inhibit_x.get(0);
            for (int i = 1; i < inhibit_x.size(); i++) {
                result += "," + inhibit_x.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write("> Activate transporter genes from " + drug_x + "\n");
        activate_x = (LinkedList) activate_table.get(drug_x);
        if (activate_x.size() == 0) {
            result = "0";
        } else {
            result = activate_x.size() + "\t" + activate_x.get(0);
            for (int i = 1; i < activate_x.size(); i++) {
                result += "," + activate_x.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write("> Inhibit transporter genes from " + drug_y + "\n");
        inhibit_y = (LinkedList) inhibit_table.get(drug_y);
        if (inhibit_y.size() == 0) {
            result = "0";
        } else {
            result = inhibit_y.size() + "\t" + inhibit_y.get(0);
            for (int i = 1; i < inhibit_y.size(); i++) {
                result += "," + inhibit_y.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write("> Activate transporter genes from " + drug_y + "\n");
        activate_y = (LinkedList) activate_table.get(drug_y);
        if (activate_y.size() == 0) {
            result = "0";
        } else {
            result = activate_y.size() + "\t" + activate_y.get(0);
            for (int i = 1; i < activate_y.size(); i++) {
                result += "," + activate_y.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write("> Common genes from both drugs\n");
        bw.write(">> same direction, sensitive genes\n");
        if (common_same_d_sensitive.size() == 0) {
            result = "0";
        } else {
            result = common_same_d_sensitive.size() + "\t" + common_same_d_sensitive.get(0);
            for (int i = 1; i < common_same_d_sensitive.size(); i++) {
                result += "," + common_same_d_sensitive.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write(">> same direction, non_sensitive genes\n");
        if (common_same_d_non_sensitive.size() == 0) {
            result = "0";
        } else {
            result = common_same_d_non_sensitive.size() + "\t" + common_same_d_non_sensitive.get(0);
            for (int i = 1; i < common_same_d_non_sensitive.size(); i++) {
                result += "," + common_same_d_non_sensitive.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write(">> opposite direction, sensitive genes\n");
        if (common_opposite_d_sensitive.size() == 0) {
            result = "0";
        } else {
            result = common_opposite_d_sensitive.size() + "\t" + common_opposite_d_sensitive.get(0);
            for (int i = 1; i < common_opposite_d_sensitive.size(); i++) {
                result += "," + common_opposite_d_sensitive.get(i);
            }
        }
        bw.write(result + "\n");

        bw.write(">> opposite direction, non_sensitive genes\n");
        if (common_opposite_d_non_sensitive.size() == 0) {
            result = "0";
        } else {
            result = common_opposite_d_non_sensitive.size() + "\t" + common_opposite_d_non_sensitive.get(0);
            for (int i = 1; i < common_opposite_d_non_sensitive.size(); i++) {
                result += "," + common_opposite_d_non_sensitive.get(i);
            }
        }
        bw.write(result + "\n");

        bw.flush();
        bw.close();
        fw.close();
    }

    public void makeDrugCombinationScore_with_tr_gene(String drug_name_list, String input_dir, String tr_gene_input_dir, String output) throws FileNotFoundException, IOException {   // 121001
        FileReader fr = new FileReader(drug_name_list);
        BufferedReader br = new BufferedReader(fr);

        FileWriter fw = new FileWriter(output);
        BufferedWriter bw = new BufferedWriter(fw);
        bw.write("combination\tscore\tsynergy\tantagonist\tdenominator\n");

        LinkedList drug_list = new LinkedList();

        String temp = br.readLine();
        while (temp != null) {
            drug_list.add(temp);
            temp = br.readLine();
        }
        br.close();
        fr.close();
        Collections.sort(drug_list);

        for (int i = 0; i < drug_list.size() - 1; i++) {
            for (int j = i + 1; j < drug_list.size(); j++) {
                String combination = (String) drug_list.get(i) + "_" + (String) drug_list.get(j);
                String filename = input_dir + "combination/" + combination + "_combination_stats.txt";
                DrugCombination dc = new DrugCombination();
                dc.makeDrugCombinationProfile(filename);
                double score = dc.getCombinationScore_with_tr_gene(input_dir, tr_gene_input_dir);
//                double score = dc.getCombinationScore_with_tr_gene_sqrt(input_dir, tr_gene_input_dir);
//                double score = dc.getCombinationScore_with_tr_gene_reverse(input_dir, tr_gene_input_dir);
//                double score = dc.getCombinationScore_full(input_dir, tr_gene_input_dir);
//                double score = dc.getCombinationScore_with_t_diff(input_dir, tr_gene_input_dir);
                double synergy = dc.getSynergy_by_tscore(input_dir);
                double antagonist = dc.getAntagonist_term(tr_gene_input_dir);
                double denominator = dc.getRevisedDenominator();

                bw.write(combination + "\t" + score + "\t" + synergy + "\t" + antagonist + "\t" + denominator + "\n");
            }
        }

        bw.flush();
        bw.close();
        fw.close();
    }

    public void makeDrugCombinationScore_with_tr_gene_antagonist_separate(String drug_name_list, String input_dir, String tr_gene_input_dir, String output) throws FileNotFoundException, IOException {   // 121001
        FileReader fr = new FileReader(drug_name_list);
        BufferedReader br = new BufferedReader(fr);

        FileWriter fw = new FileWriter(output);
        BufferedWriter bw = new BufferedWriter(fw);
        bw.write("combination\tscore\tsynergy\tantagonist\tdenominator\n");

        LinkedList drug_list = new LinkedList();

        String temp = br.readLine();
        while (temp != null) {
            drug_list.add(temp);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        for (int i = 0; i < drug_list.size() - 1; i++) {
            for (int j = i + 1; j < drug_list.size(); j++) {
                String combination = (String) drug_list.get(i) + "_" + (String) drug_list.get(j);
                String filename = input_dir + "combination/" + combination + "_combination_stats.txt";
                DrugCombination dc = new DrugCombination();
                dc.makeDrugCombinationProfile(filename);
                double score = dc.getCombinationScore_with_tr_gene_antagonist_separate(input_dir, tr_gene_input_dir);
                double synergy = dc.getSynergy_by_tscore(input_dir);
                double antagonist = dc.getAntagonist_term_separate(tr_gene_input_dir);
                double denominator = dc.getRevisedDenominator();

                bw.write(combination + "\t" + score + "\t" + synergy + "\t" + antagonist + "\t" + denominator + "\n");
            }
        }

        bw.flush();
        bw.close();
        fw.close();
    }

    public void makeDrugCombinationScore_with_tr_gene_between_same_gene(String drug_name_list, String input_dir, String tr_gene_input_dir, String output) throws FileNotFoundException, IOException {   // 121001
        FileReader fr = new FileReader(drug_name_list);
        BufferedReader br = new BufferedReader(fr);

        FileWriter fw = new FileWriter(output);
        BufferedWriter bw = new BufferedWriter(fw);
        bw.write("combination\tscore\tsynergy\tantagonist\tdenominator\n");

        LinkedList drug_list = new LinkedList();

        String temp = br.readLine();
        while (temp != null) {
            drug_list.add(temp);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        for (int i = 0; i < drug_list.size(); i++) {
            String combination = (String) drug_list.get(i) + "_" + (String) drug_list.get(i);
//            String filename = input_dir + "combination/" + combination + "_combination_stats.txt";
            String filename = input_dir + "same/" + combination + "_combination_stats.txt";
            DrugCombination dc = new DrugCombination();
            dc.makeDrugCombinationProfile(filename);
            double score = dc.getCombinationScore_with_tr_gene_antagonist_separate(input_dir, tr_gene_input_dir);
            double synergy = dc.getSynergy_by_tscore(input_dir);
            double antagonist = dc.getAntagonist_term_separate(tr_gene_input_dir);
            double denominator = dc.getRevisedDenominator();

            bw.write(combination + "\t" + score + "\t" + synergy + "\t" + antagonist + "\t" + denominator + "\n");

        }

        bw.flush();
        bw.close();
        fw.close();
    }
}
