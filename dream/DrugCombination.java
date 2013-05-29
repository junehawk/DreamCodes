package dream;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class DrugCombination {

    private String name_a;
    private String name_b;
    private int only_a;
    private int only_b;
    private int same_d;
    private int opposite_d;
    private String only_a_gene;
    private String only_b_gene;
    private String same_d_gene;
    private String opposite_d_gene;
    private String inhibit_a_gene;
    private String activate_a_gene;
    private String inhibit_b_gene;
    private String activate_b_gene;
    private int inhibit_a;
    private int activate_a;
    private int inhibit_b;
    private int activate_b;

    /**
     * @return the only_a
     */
    public void makeDrugCombinationProfile(String filename) throws FileNotFoundException, IOException {
        FileReader fr = new FileReader(filename);
        BufferedReader br = new BufferedReader(fr);
        String temp = br.readLine();    // only_a_header
        String[] s;

        setName_a(temp.replace("> Genes only from ", ""));
        temp = br.readLine();
        s = temp.split("\t");
        setOnly_a(Integer.parseInt(s[0]));
        if (getOnly_a() > 0) {
            setOnly_a_gene(s[1]);
        }
        temp = br.readLine();

        setName_b(temp.replace("> Genes only from ", ""));
        temp = br.readLine();
        s = temp.split("\t");
        setOnly_b(Integer.parseInt(s[0]));
        if (getOnly_b() > 0) {
            setOnly_b_gene(s[1]);
        }
        temp = br.readLine();

        temp = br.readLine();   // inhibit for a
        s = temp.split("\t");
        setInhibit_a(Integer.parseInt(s[0]));
        if (getInhibit_a() > 0) {
            setInhibit_a_gene(s[1]);
        }
        temp = br.readLine();

        temp = br.readLine();   // activate for a
        s = temp.split("\t");
        setActivate_a(Integer.parseInt(s[0]));
        if (getActivate_a() > 0) {
            setActivate_a_gene(s[1]);
        }
        temp = br.readLine();

        temp = br.readLine();   // inhibit for b
        s = temp.split("\t");
        setInhibit_b(Integer.parseInt(s[0]));
        if (getInhibit_b() > 0) {
            setInhibit_b_gene(s[1]);
        }
        temp = br.readLine();

        temp = br.readLine();   // activate for b
        s = temp.split("\t");
        setActivate_b(Integer.parseInt(s[0]));
        if (getActivate_b() > 0) {
            setActivate_b_gene(s[1]);
        }
        temp = br.readLine();

        setSame_d(0);
        setOpposite_d(0);
        setSame_d_gene("");
        setOpposite_d_gene("");

        while (temp != null) {
            if (temp.indexOf("same") != -1) {
                temp = br.readLine();
                s = temp.split("\t");
                int cnt = Integer.parseInt(s[0]);
                if (cnt > 0) {
                    setSame_d(getSame_d() + cnt);
                    if (getSame_d_gene().equals("")) {
                        setSame_d_gene(s[1]);
                    } else {
                        setSame_d_gene(getSame_d_gene() + "," + s[1]);
                    }
                }
            }
            if (temp.indexOf("opposite") != -1) {
                temp = br.readLine();
                s = temp.split("\t");
                int cnt = Integer.parseInt(s[0]);
                if (cnt > 0) {
                    setOpposite_d(getOpposite_d() + cnt);
                    if (getOpposite_d_gene().equals("")) {
                        setOpposite_d_gene(s[1]);
                    } else {
                        setOpposite_d_gene(getOpposite_d_gene() + "," + s[1]);
                    }
                }
            }
            temp = br.readLine();
        }
    }

    public double getCombinationScore(String input_dir) throws FileNotFoundException, IOException {
        double score;
        double denominator;
        double total_a = getTotal_a();
        double total_b = getTotal_b();
        if (total_a != 0 && total_b != 0) {
            denominator = ((double) getOnly_a() / total_a) + ((double) getOnly_b() / total_b);
        } else if (total_a == 0) {
            denominator = 1 + ((double) getOnly_b() / total_b);
        } else {
            denominator = ((double) getOnly_a() / total_a) + 1;
        }

        double numerator = getSynergy_by_tscore(input_dir);
        score = numerator / denominator;

        return score;
    }

    public double getSynergy_by_tscore(String input_dir) throws FileNotFoundException, IOException {
        double score = 0;
        String[] s;
        Hashtable a_table = new Hashtable(10000);
        Hashtable b_table = new Hashtable(10000);

        String filename = input_dir + getName_a() + "_24_IC20_gene_based_sig_filtered.txt";
        FileReader fr = new FileReader(filename);
        BufferedReader br = new BufferedReader(fr);
        String temp = br.readLine();    // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            a_table.put(s[0], s[1]);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        filename = input_dir + getName_b() + "_24_IC20_gene_based_sig_filtered.txt";
        fr = new FileReader(filename);
        br = new BufferedReader(fr);
        temp = br.readLine();   // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            b_table.put(s[0], s[1]);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        s = getSame_d_gene().split(",");
        if (s.length > 1) {
            for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                System.out.println(s[i]);
//                System.out.println(getSame_d_gene());
                double score_a = Double.parseDouble((String) a_table.get(s[i]));
                double score_b = Double.parseDouble((String) b_table.get(s[i]));
                score += score_a * score_b;
            }
        }

        s = getOpposite_d_gene().split(",");
        if (s.length > 1) {
            for (int i = 0; i < s.length; i++) {
                double score_a = Double.parseDouble((String) a_table.get(s[i]));
                double score_b = Double.parseDouble((String) b_table.get(s[i]));
                score += score_a * score_b;
            }
        }

        if (getSame_d() + getOpposite_d() == 0) {
            score = 0;
        } else {
            score = score / (getSame_d() + getOpposite_d());
        }

        return score;
    }

    public double getSynergy_by_tscore_diff(String input_dir) throws FileNotFoundException, IOException {
        double score = 0;
        String[] s;
        Hashtable a_table = new Hashtable(10000);
        Hashtable b_table = new Hashtable(10000);

        String filename = input_dir + getName_a() + "_24_IC20_gene_based_sig_filtered.txt";
        FileReader fr = new FileReader(filename);
        BufferedReader br = new BufferedReader(fr);
        String temp = br.readLine();    // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            a_table.put(s[0], s[1]);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        filename = input_dir + getName_b() + "_24_IC20_gene_based_sig_filtered.txt";
        fr = new FileReader(filename);
        br = new BufferedReader(fr);
        temp = br.readLine();   // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            b_table.put(s[0], s[1]);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        s = getSame_d_gene().split(",");
        if (s.length > 1) {
            for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                System.out.println(s[i]);
//                System.out.println(getSame_d_gene());
//                double score_a = Double.parseDouble((String) a_table.get(s[i]));
//                double score_b = Double.parseDouble((String) b_table.get(s[i]));
//                score += score_a * score_b;
            }
        }

        s = getOpposite_d_gene().split(",");
        if (s.length > 1) {
            for (int i = 0; i < s.length; i++) {
                double score_a = Double.parseDouble((String) a_table.get(s[i]));
                double score_b = Double.parseDouble((String) b_table.get(s[i]));
                score += score_a * score_b;
            }
        }

        if (getSame_d() + getOpposite_d() == 0) {
            score = 0;
        } else {
            score = score / getOpposite_d();
        }

        return score;
    }

    public double getSynergy_by_tscore_abs(String input_dir) throws FileNotFoundException, IOException {
        double score = 0;
        String[] s;
        Hashtable a_table = new Hashtable(10000);
        Hashtable b_table = new Hashtable(10000);

        String filename = input_dir + getName_a() + "_24_IC20_gene_based_sig_filtered.txt";
        FileReader fr = new FileReader(filename);
        BufferedReader br = new BufferedReader(fr);
        String temp = br.readLine();    // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            a_table.put(s[0], s[1]);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        filename = input_dir + getName_b() + "_24_IC20_gene_based_sig_filtered.txt";
        fr = new FileReader(filename);
        br = new BufferedReader(fr);
        temp = br.readLine();   // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            b_table.put(s[0], s[1]);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        s = getSame_d_gene().split(",");
        if (s.length > 1) {
            for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                System.out.println(s[i]);
//                System.out.println(getSame_d_gene());
                double score_a = Double.parseDouble((String) a_table.get(s[i]));
                double score_b = Double.parseDouble((String) b_table.get(s[i]));
                score += Math.abs(score_a * score_b);
            }
        }

        s = getOpposite_d_gene().split(",");
        if (s.length > 1) {
            for (int i = 0; i < s.length; i++) {
                double score_a = Double.parseDouble((String) a_table.get(s[i]));
                double score_b = Double.parseDouble((String) b_table.get(s[i]));
                score += Math.abs(score_a * score_b);
            }
        }

        if (getSame_d() + getOpposite_d() == 0) {
            score = 0;
        } else {
            score = score / (getSame_d() + getOpposite_d());
        }

        return score;
    }

    public double getSynergy_by_only_gene(String input_dir) throws FileNotFoundException, IOException {
        double score_a = 0;
        double score_b = 0;
        String[] s;
        Hashtable a_table = new Hashtable(10000);
        Hashtable b_table = new Hashtable(10000);

        String filename = input_dir + getName_a() + "_24_IC20_gene_based_sig_filtered.txt";
        FileReader fr = new FileReader(filename);
        BufferedReader br = new BufferedReader(fr);
        String temp = br.readLine();    // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            a_table.put(s[0], s[1]);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        filename = input_dir + getName_b() + "_24_IC20_gene_based_sig_filtered.txt";
        fr = new FileReader(filename);
        br = new BufferedReader(fr);
        temp = br.readLine();   // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            b_table.put(s[0], s[1]);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        if (only_a > 0) {
            s = getOnly_a_gene().split(",");
            if (s.length > 1) {
                for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                System.out.println(s[i]);
//                System.out.println(getSame_d_gene());
                    double t_score_a = Double.parseDouble((String) a_table.get(s[i]));
                    score_a += t_score_a * t_score_a;
                }
            }
        }

        if (only_b > 0) {
            s = getOnly_b_gene().split(",");
            if (s.length > 1) {
                for (int i = 0; i < s.length; i++) {
                    double t_score_b = Double.parseDouble((String) b_table.get(s[i]));
                    score_b += t_score_b * t_score_b;
                }
            }
            if (only_a != 0) {
                score_a = score_a / only_a;
            }
            if (only_b != 0) {
                score_b = score_b / only_b;
            }
        }

        return score_a + score_b;
    }

    public double getCombinationScore_with_tr_gene(String input_dir, String tr_gene_input_dir) throws FileNotFoundException, IOException {
        double score;
        double denominator = getRevisedDenominator();
        double synergy_term = getSynergy_by_tscore(input_dir);
        double antagonist_term = getAntagonist_term(tr_gene_input_dir);

        score = (synergy_term + antagonist_term) / denominator;

        return score;
    }

    public double getCombinationScore_full(String input_dir, String tr_gene_input_dir) throws FileNotFoundException, IOException {
        double score;
        double denominator = getRevisedDenominator();
        double only_term = getSynergy_by_only_gene(input_dir);
        double synergy_term = getSynergy_by_tscore(input_dir);
        double antagonist_term = getAntagonist_term(tr_gene_input_dir);

        score = (only_term - synergy_term + antagonist_term) / denominator;

        return score;
    }

    public double getCombinationScore_with_t_diff(String input_dir, String tr_gene_input_dir) throws FileNotFoundException, IOException {
        double score;
        double denominator = getRevisedDenominator();
        double synergy_term = getSynergy_by_tscore_diff(input_dir);
        double antagonist_term = getAntagonist_term(tr_gene_input_dir);

        score = (synergy_term + antagonist_term) / denominator;

        return score;
    }

    public double getCombinationScore_with_tr_gene_reverse(String input_dir, String tr_gene_input_dir) throws FileNotFoundException, IOException {
        double score;
        double denominator = getRevisedDenominator();
        double synergy_term = getSynergy_by_tscore(input_dir);
//        double synergy_term = getSynergy_by_tscore_abs(input_dir);
        double antagonist_term = getAntagonist_term(tr_gene_input_dir);

        score = (-synergy_term + antagonist_term) / denominator;

        return score;
    }

    public double getCombinationScore_with_tr_gene_sqrt(String input_dir, String tr_gene_input_dir) throws FileNotFoundException, IOException {
        double score;
        double denominator = getRevisedDenominator();
        double synergy_term = getSynergy_by_tscore(input_dir);
        if (synergy_term < 0) {
            synergy_term = -1 * synergy_term;
            synergy_term = Math.sqrt(synergy_term);
            synergy_term = -1 * synergy_term;
        } else {
            synergy_term = Math.sqrt(synergy_term);
        }
        double antagonist_term = getAntagonist_term(tr_gene_input_dir);
        if (antagonist_term < 0) {
            antagonist_term = -1 * antagonist_term;
            antagonist_term = Math.sqrt(antagonist_term);
            antagonist_term = -1 * antagonist_term;
        } else {
            antagonist_term = Math.sqrt(antagonist_term);
        }

        score = (synergy_term + antagonist_term) / denominator;

        return score;
    }

    public double getCombinationScore_with_tr_gene_antagonist_separate(String input_dir, String tr_gene_input_dir) throws FileNotFoundException, IOException {
        double score;
        double denominator = getRevisedDenominator();
        double synergy_term = getSynergy_by_tscore(input_dir);
        double antagonist_term = getAntagonist_term_separate(tr_gene_input_dir);

        score = (synergy_term + antagonist_term) / denominator;

        return score;
    }

    public double getRevisedDenominator() {
        double denominator;
        double total_a = getTotal_a();
        double total_b = getTotal_b();
        if (total_a != 0 && total_b != 0) {
            denominator = ((double) (total_a - (double) (same_d + opposite_d) / 2) / total_a) + ((double) (total_b - (double) (same_d + opposite_d) / 2) / total_b);
        } else if (total_a == 0) {
            denominator = 1 + ((double) (total_b - (double) (same_d + opposite_d) / 2) / total_b);
        } else {
            denominator = ((double) (total_a - (double) (same_d + opposite_d) / 2) / total_a) + 1;
        }
        return denominator;
    }

    public double getAntagonist_term(String tr_gene_input_dir) throws FileNotFoundException, IOException {
        double score = 0;
        String[] s;
        Hashtable a_table = new Hashtable(10000);
        Hashtable b_table = new Hashtable(10000);

        String filename = tr_gene_input_dir + getName_a() + "_24_IC20_gene_based_tr_gene_filtered.txt";
        FileReader fr = new FileReader(filename);
        BufferedReader br = new BufferedReader(fr);
        String temp = br.readLine();    // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            a_table.put(s[0], temp);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        filename = tr_gene_input_dir + getName_b() + "_24_IC20_gene_based_tr_gene_filtered.txt";
        fr = new FileReader(filename);
        br = new BufferedReader(fr);
        temp = br.readLine();   // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            b_table.put(s[0], temp);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        double tr_a_score = 0;
        int num_tr_a = 0;
        if (getInhibit_a() > 0) {
            s = getInhibit_a_gene().split(",");
            if (s.length > 1) {
                for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                System.out.println(s[i]);
                    String b_stat = (String) b_table.get(s[i]);
                    if (b_stat != null) {
                        String[] g = b_stat.split("\t");
                        double flag = 1;
                        if (Double.parseDouble(g[4]) < 0.1) { // DEG
                            double t = Double.parseDouble(g[1]);
                            if (t > 0) {
                                flag = -1;  // inhibitory
                            }
                            t = flag * t * t;
                            tr_a_score += t;
                            num_tr_a++;
                        }
                    }
                }
            }
        }
        if (getActivate_a() > 0) {
            s = getActivate_a_gene().split(",");
            if (s.length > 1) {
                for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                System.out.println(s[i]);
                    String b_stat = (String) b_table.get(s[i]);
                    if (b_stat != null) {
                        String[] g = b_stat.split("\t");
                        double flag = 1;
                        if (Double.parseDouble(g[4]) < 0.1) { // DEG
                            double t = Double.parseDouble(g[1]);
                            if (t < 0) {
                                flag = -1;  // activative
                            }
                            t = flag * t * t;
                            tr_a_score += t;
                            num_tr_a++;
                        }
                    }
                }
            }
        }

        double tr_b_score = 0;
        int num_tr_b = 0;
        if (getInhibit_b() > 0) {
            s = getInhibit_b_gene().split(",");
            if (s.length > 1) {
                for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                    System.out.println(s[i]);
                    String a_stat = (String) a_table.get(s[i]);
                    if (a_stat != null) {
                        String[] g = a_stat.split("\t");
                        double flag = 1;
                        if (Double.parseDouble(g[4]) < 0.1) { // DEG
                            double t = Double.parseDouble(g[1]);
                            if (t > 0) {
                                flag = -1;  // inhibitory
                            }
                            t = flag * t * t;
                            tr_b_score += t;
                            num_tr_b++;
                        }
                    }
                }
            }
        }
        if (getActivate_b() > 0) {
            s = getActivate_b_gene().split(",");
            if (s.length > 1) {
                for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                System.out.println(s[i]);
                    String a_stat = (String) a_table.get(s[i]);
                    if (a_stat != null) {
                        String[] g = a_stat.split("\t");
                        double flag = 1;
                        if (Double.parseDouble(g[4]) < 0.1) { // DEG
                            double t = Double.parseDouble(g[1]);
                            if (t < 0) {
                                flag = -1;  // activative
                            }
                            t = flag * t * t;
                            tr_b_score += t;
                            num_tr_b++;
                        }
                    }
                }
            }
        }

        if (num_tr_a + num_tr_b != 0) {
            score = (double) (tr_a_score + tr_b_score) / (num_tr_a + num_tr_b);
        }

        return score;
    }

    public double getAntagonist_term_separate(String tr_gene_input_dir) throws FileNotFoundException, IOException {
        double score = 0;
        String[] s;
        Hashtable a_table = new Hashtable(10000);
        Hashtable b_table = new Hashtable(10000);

        String filename = tr_gene_input_dir + getName_a() + "_24_IC20_gene_based_tr_gene_filtered.txt";
        FileReader fr = new FileReader(filename);
        BufferedReader br = new BufferedReader(fr);
        String temp = br.readLine();    // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            a_table.put(s[0], temp);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        filename = tr_gene_input_dir + getName_b() + "_24_IC20_gene_based_tr_gene_filtered.txt";
        fr = new FileReader(filename);
        br = new BufferedReader(fr);
        temp = br.readLine();   // header
        temp = br.readLine();
        while (temp != null) {
            s = temp.split("\t");
            b_table.put(s[0], temp);
            temp = br.readLine();
        }
        br.close();
        fr.close();

        double tr_a_score = 0;
        int num_tr_a = 0;
        if (getInhibit_a() > 0) {
            s = getInhibit_a_gene().split(",");
            if (s.length > 1) {
                for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                System.out.println(s[i]);
                    String b_stat = (String) b_table.get(s[i]);
                    if (b_stat != null) {
                        String[] g = b_stat.split("\t");
                        double flag = 1;
                        if (Double.parseDouble(g[4]) < 0.1) { // DEG
                            double t = Double.parseDouble(g[1]);
                            if (t > 0) {
                                flag = -1;  // inhibitory
                            }
                            t = flag * t * t;
                            tr_a_score += t;
                            num_tr_a++;
                        }
                    }
                }
            }
        }
        if (getActivate_a() > 0) {
            s = getActivate_a_gene().split(",");
            if (s.length > 1) {
                for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                System.out.println(s[i]);
                    String b_stat = (String) b_table.get(s[i]);
                    if (b_stat != null) {
                        String[] g = b_stat.split("\t");
                        double flag = 1;
                        if (Double.parseDouble(g[4]) < 0.1) { // DEG
                            double t = Double.parseDouble(g[1]);
                            if (t < 0) {
                                flag = -1;  // activative
                            }
                            t = flag * t * t;
                            tr_a_score += t;
                            num_tr_a++;
                        }
                    }
                }
            }
        }

        double tr_b_score = 0;
        int num_tr_b = 0;
        if (getInhibit_b() > 0) {
            s = getInhibit_b_gene().split(",");
            if (s.length > 1) {
                for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                    System.out.println(s[i]);
                    String a_stat = (String) a_table.get(s[i]);
                    if (a_stat != null) {
                        String[] g = a_stat.split("\t");
                        double flag = 1;
                        if (Double.parseDouble(g[4]) < 0.1) { // DEG
                            double t = Double.parseDouble(g[1]);
                            if (t > 0) {
                                flag = -1;  // inhibitory
                            }
                            t = flag * t * t;
                            tr_b_score += t;
                            num_tr_b++;
                        }
                    }
                }
            }
        }
        if (getActivate_b() > 0) {
            s = getActivate_b_gene().split(",");
            if (s.length > 1) {
                for (int i = 0; i < s.length; i++) {
//                System.out.println(name_a+"\t"+name_b);
//                System.out.println(s[i]);
                    String a_stat = (String) a_table.get(s[i]);
                    if (a_stat != null) {
                        String[] g = a_stat.split("\t");
                        double flag = 1;
                        if (Double.parseDouble(g[4]) < 0.1) { // DEG
                            double t = Double.parseDouble(g[1]);
                            if (t < 0) {
                                flag = -1;  // activative
                            }
                            t = flag * t * t;
                            tr_b_score += t;
                            num_tr_b++;
                        }
                    }
                }
            }
        }

        if (num_tr_a != 0 && num_tr_b != 0) {
            score = (double) tr_a_score / num_tr_a + (double) tr_b_score / num_tr_b;
        } else if (num_tr_a == 0 && num_tr_b != 0) {
            score = (double) tr_b_score / num_tr_b;
        } else if (num_tr_a != 0 && num_tr_b == 0) {
            score = (double) tr_a_score / num_tr_a;
        }

        return score;
    }

    public int getOnly_a() {
        return only_a;
    }

    /**
     * @param only_a the only_a to set
     */
    public void setOnly_a(int only_a) {
        this.only_a = only_a;
    }

    /**
     * @return the only_b
     */
    public int getOnly_b() {
        return only_b;
    }

    /**
     * @param only_b the only_b to set
     */
    public void setOnly_b(int only_b) {
        this.only_b = only_b;
    }

    /**
     * @return the same_d
     */
    public int getSame_d() {
        return same_d;
    }

    /**
     * @param same_d the same_d to set
     */
    public void setSame_d(int same_d) {
        this.same_d = same_d;
    }

    /**
     * @return the opposite_d
     */
    public int getOpposite_d() {
        return opposite_d;
    }

    /**
     * @param opposite_d the opposite_d to set
     */
    public void setOpposite_d(int opposite_d) {
        this.opposite_d = opposite_d;
    }

    /**
     * @return the only_a_gene
     */
    public String getOnly_a_gene() {
        return only_a_gene;
    }

    /**
     * @param only_a_gene the only_a_gene to set
     */
    public void setOnly_a_gene(String only_a_gene) {
        this.only_a_gene = only_a_gene;
    }

    /**
     * @return the only_b_gene
     */
    public String getOnly_b_gene() {
        return only_b_gene;
    }

    /**
     * @param only_b_gene the only_b_gene to set
     */
    public void setOnly_b_gene(String only_b_gene) {
        this.only_b_gene = only_b_gene;
    }

    /**
     * @return the same_d_gene
     */
    public String getSame_d_gene() {
        return same_d_gene;
    }

    /**
     * @param same_d_gene the same_d_gene to set
     */
    public void setSame_d_gene(String same_d_gene) {
        this.same_d_gene = same_d_gene;
    }

    /**
     * @return the opposite_d_gene
     */
    public String getOpposite_d_gene() {
        return opposite_d_gene;
    }

    /**
     * @param opposite_d_gene the opposite_d_gene to set
     */
    public void setOpposite_d_gene(String opposite_d_gene) {
        this.opposite_d_gene = opposite_d_gene;
    }

    /**
     * @return the name_a
     */
    public String getName_a() {
        return name_a;
    }

    /**
     * @param name_a the name_a to set
     */
    public void setName_a(String name_a) {
        this.name_a = name_a;
    }

    /**
     * @return the name_b
     */
    public String getName_b() {
        return name_b;
    }

    /**
     * @param name_b the name_b to set
     */
    public void setName_b(String name_b) {
        this.name_b = name_b;
    }

    public int getTotal_a() {
        return getOnly_a() + getSame_d() + getOpposite_d();
    }

    public int getTotal_b() {
        return getOnly_b() + getSame_d() + getOpposite_d();
    }

    /**
     * @return the inhibit_a_gene
     */
    public String getInhibit_a_gene() {
        return inhibit_a_gene;
    }

    /**
     * @param inhibit_a_gene the inhibit_a_gene to set
     */
    public void setInhibit_a_gene(String inhibit_a_gene) {
        this.inhibit_a_gene = inhibit_a_gene;
    }

    /**
     * @return the activate_a_gene
     */
    public String getActivate_a_gene() {
        return activate_a_gene;
    }

    /**
     * @param activate_a_gene the activate_a_gene to set
     */
    public void setActivate_a_gene(String activate_a_gene) {
        this.activate_a_gene = activate_a_gene;
    }

    /**
     * @return the inhibit_b_gene
     */
    public String getInhibit_b_gene() {
        return inhibit_b_gene;
    }

    /**
     * @param inhibit_b_gene the inhibit_b_gene to set
     */
    public void setInhibit_b_gene(String inhibit_b_gene) {
        this.inhibit_b_gene = inhibit_b_gene;
    }

    /**
     * @return the activate_b_gene
     */
    public String getActivate_b_gene() {
        return activate_b_gene;
    }

    /**
     * @param activate_b_gene the activate_b_gene to set
     */
    public void setActivate_b_gene(String activate_b_gene) {
        this.activate_b_gene = activate_b_gene;
    }

    /**
     * @return the inhibit_a
     */
    public int getInhibit_a() {
        return inhibit_a;
    }

    /**
     * @param inhibit_a the inhibit_a to set
     */
    public void setInhibit_a(int inhibit_a) {
        this.inhibit_a = inhibit_a;
    }

    /**
     * @return the activate_a
     */
    public int getActivate_a() {
        return activate_a;
    }

    /**
     * @param activate_a the activate_a to set
     */
    public void setActivate_a(int activate_a) {
        this.activate_a = activate_a;
    }

    /**
     * @return the inhibit_b
     */
    public int getInhibit_b() {
        return inhibit_b;
    }

    /**
     * @param inhibit_b the inhibit_b to set
     */
    public void setInhibit_b(int inhibit_b) {
        this.inhibit_b = inhibit_b;
    }

    /**
     * @return the activate_b
     */
    public int getActivate_b() {
        return activate_b;
    }

    /**
     * @param activate_b the activate_b to set
     */
    public void setActivate_b(int activate_b) {
        this.activate_b = activate_b;
    }
}
