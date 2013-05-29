/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package dream;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 *
 * @author kimjh
 */
public class Mutation {
    private Hashtable gene_table;
    
    public Mutation(String gene_list) throws FileNotFoundException, IOException{
        gene_table = new Hashtable(20000);
        FileReader fr = new FileReader(gene_list);
        BufferedReader br = new BufferedReader(fr);
        String temp = br.readLine();
        while(temp!=null){
            if(!temp.equals("")){
                gene_table.put(temp,"0");
            }
            temp = br.readLine();
        }
    }

    /**
     * @return the gene_table
     */
    public Hashtable getGene_table() {
        return gene_table;
    }

    /**
     * @param gene_table the gene_table to set
     */
    public void setGene_table(Hashtable gene_table) {
        this.gene_table = gene_table;
    }
    
    public void setMutation(String gene){
        this.gene_table.put(gene,"1");
    }
    
    public LinkedList getSortedGeneList(){
        Set set = gene_table.keySet();
        LinkedList gene_list = new LinkedList(set);
        Collections.sort(gene_list);
        return gene_list;
    }
    
    public LinkedList getSortedMutationStatus(){
        Set set = gene_table.keySet();
        LinkedList gene_list = new LinkedList(set);
        Collections.sort(gene_list);
        Iterator keys = gene_list.iterator();
        LinkedList mutation_status = new LinkedList();
        while(keys.hasNext()){
            String key = (String) keys.next();
            String mut = (String) gene_table.get(key);
            mutation_status.add(mut);
        }
        return mutation_status;
    }
    
    public int getNumberOfGene(){
        return gene_table.size();
    }
}
