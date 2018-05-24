/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utils;
import contiguousfspm.pseudoSequentialPattern;
import htsjdk.samtools.QueryInterval;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author jiadonglin
 */
public class svOutInfo {
     
    int start;
    int end;
    String pattern;
    int linkType; // 1 mate linked. -1 self linked. 0 non-linkable
    int[] supEvi;
    int[] selfLinkedMapQ;
    int[] selfLinkedWeight;
    String[] selfLinkedItemTypes;
    // Used for saving SVs called from a pattern, usually > 100bp
    List<Integer> weights;
    List<Integer> postions;    
    List<Double> ratios;
    List<String> oris;
    int[] arpSpanBpMapQ;
    int[] arpSpanBpWeight;
    String[] arpSpanBpItem;
    
    boolean isPassed = true;
    
    
    
    public svOutInfo(int s, int e, String patternStr, int linkFlag, int[] supEvidence, List<Integer> weight, List<Integer> pos, List<Double> wRatio, List<String> oris){        
        start = s;
        end = e;

        if (s > e && e != -1){
            start = e;
            end = s;
        }
        pattern = patternStr;
        linkType = linkFlag;
        supEvi = supEvidence;
        weights = weight;
        postions = pos;
        ratios = wRatio;
        this.oris = oris;
    }
    
    public svOutInfo(int s, int e, int linkFlag, int weight, int plusReadNum, int minusReadNum){
        start = s;
        end = e;
        supEvi = new int[3];
        supEvi[0] = weight;
        supEvi[1] = plusReadNum;
        supEvi[2] = minusReadNum;
        linkType = linkFlag;

    }
    @Override
    public String toString(){
        if (linkType == 9 || linkType == -3 || linkType == -4 || linkType == -5){
            isPassed = false;
        }
        StringBuilder sb = new StringBuilder();
        sb.append("\t");
        if (start == -1 && end != -1){
            start = end;
            end = -1;
        }
        sb.append(start);
        sb.append("\t");
        if (end == -1){
            sb.append("-");
        }else{
            sb.append(end);
        }
                
        sb.append("\t");
        if (isPassed){
            sb.append("PASS");
        }else{
            sb.append("LowQual");
        }
        sb.append("\t");
        sb.append(infoToString());
        return sb.toString();
    }
    private String infoToString(){
        StringBuilder sb = new StringBuilder();
        sb.append("SupType=");
        if (linkType == -5){
            sb.append("MultiBP");
        }
        if (linkType == -4){
            sb.append("OEM;");
            sb.append("OEMWeight=");
            sb.append(supEvi[1]);
        }
        
        if (linkType == -3){
            sb.append("None");
        }
        
        if (linkType == -2){
            sb.append("SMALL_INSERT");
        }
        
        if(linkType == -1){
            sb.append("Self");                        
            sb.append(";LinkedMapQ=");
            sb.append(intArrayToString(selfLinkedMapQ));
            sb.append(";LinkedWeight=");
            sb.append(intArrayToString(selfLinkedWeight));
            sb.append(";LinkedItem=");
            sb.append(strArrayToString(selfLinkedItemTypes));
        }
        if(linkType == 0){
            sb.append("None");            
            sb.append(";SUP=NA");
        }
        if (linkType == 1){
            sb.append("ARP_Span");               
            sb.append(";ARP_SUP=");
            sb.append(supEvi[0]);
            sb.append(";BP_mapQ=");
            sb.append(intArrayToString(arpSpanBpMapQ));
            sb.append(";BP_weight=");
            sb.append(intArrayToString(arpSpanBpWeight));
            sb.append(";BP_item=");
            sb.append(strArrayToString(arpSpanBpItem));
        }                
        if (linkType == 2){
            sb.append("ARP_Span&Split;");            
            sb.append("ARP_SUP=");
            sb.append(supEvi[0]);
            sb.append(";Split_SUP=");
            sb.append(supEvi[3]);
            
        }
        if (linkType == 3){
            sb.append("Split;");
            sb.append("Split_SUP=");
            sb.append(supEvi[3]);
            sb.append(";Split_mapQ=");
            sb.append(supEvi[4]);
        }
        if (linkType == 4){
            sb.append("Self&Cross;");
            if (supEvi[0] != 0){
                sb.append("ARP_SHARE=");
                sb.append(supEvi[0]);                
            }            
            sb.append(";CROSS_LEN=");
            sb.append(supEvi[1]);            
            sb.append(";CROSS_READ=");
            sb.append(supEvi[2]);
            
        }
        if (linkType == 5){
            sb.append("Self&Split;");            
            sb.append("ARP_SHARE=");
            sb.append(supEvi[0]);
            sb.append(";Split_SUP=");
            sb.append(supEvi[3]);
            sb.append(";Split_mapQ=");
            sb.append(supEvi[4]);
        }
        if (linkType == 6){
            sb.append("Self&Split&Cross;");
            if (supEvi[0] != 0){
                sb.append("ARP_SHARE=");
                sb.append(supEvi[0]);                
                sb.append(";Split_SUP=");
                sb.append(supEvi[3]);
                sb.append(";Split_mapQ=");
                sb.append(supEvi[4]);
            }             
            sb.append(";CROSS_LEN=");
            sb.append(supEvi[1]);            
            sb.append(";CROSS_READ=");
            sb.append(supEvi[2]);
        }
        if (linkType == 7){
            sb.append("Split&Cross;"); 
            sb.append("Split_SUP=");
            sb.append(supEvi[3]);
            sb.append(";Split_mapQ=");
            sb.append(supEvi[4]);
            sb.append("CROSS_LEN=");
            sb.append(supEvi[1]);            
            sb.append(";CROSS_READ=");
            sb.append(supEvi[2]);
        }
        if (linkType == 8){
            sb.append("Cross");            
            sb.append(";CROSS_LEN=");
            sb.append(supEvi[1]);            
            sb.append(";CROSS_READ=");
            sb.append(supEvi[2]);
        }
        if(linkType == 9){
            sb.append("ARP_INFER;");
            sb.append("ARP_SUP=");
            sb.append(supEvi[1]);
        }
        if (linkType == 10){
            sb.append("Realign;");
            sb.append("AR=");
            sb.append(supEvi[0]);                      
            sb.append(";PR=");
            sb.append(supEvi[1]);            
            sb.append(";MR=");
            sb.append(supEvi[2]);
        }
        
        sb.append(";Pattern=");
        sb.append(pattern);                                             
        sb.append(";Weights=");
        sb.append(intListToString(weights));        
        sb.append(";Pos=");
        sb.append(intListToString(postions));
        sb.append(";Ratio=");
        sb.append(doubleListToString(ratios));
        sb.append(";Ori=");
        sb.append(strListToString(oris));
        return sb.toString();
    }
    private String intListToString(List<Integer> alist){
        StringBuilder sb = new StringBuilder();
        for (Integer ele : alist){
            sb.append(ele);
            sb.append(",");
        }
        String str = sb.toString();
        String outStr = str.substring(0, str.length() - 1);
        return outStr;
    }
    private String doubleListToString(List<Double> alist){
        StringBuilder sb = new StringBuilder();
        for (Double ele : alist){
            sb.append(ele);
            sb.append(",");
        }
        String str = sb.toString();
        String outStr = str.substring(0, str.length() - 1);
        return outStr;
    }
    private String strListToString(List<String> alist){
        StringBuilder sb = new StringBuilder();
        for (String ele : alist){
            sb.append(ele);
            sb.append(",");
        }
        String str = sb.toString();
        String outStr = str.substring(0, str.length() - 1);
        return outStr;
    }
    
    private String intArrayToString(int[] array){
        StringBuilder sb = new StringBuilder();
        for (int ele : array){
            sb.append(ele);
            sb.append(",");
        }
        String str = sb.toString();
        String outStr = str.substring(0, str.length() - 1);
        return outStr;
    }
    
    private String strArrayToString(String[] array){
        StringBuilder sb = new StringBuilder();
        for (String ele : array){
            sb.append(ele);
            sb.append(",");
        }
        String str = sb.toString();
        String outStr = str.substring(0, str.length() - 1);
        return outStr;
    }
    
    public void setSvInfo(String patternStr, List<Integer> weights, List<Integer> posList, List<Double> ratios, List<String> oris){
        pattern = patternStr;
        this.weights = weights;
        this.postions = posList;
        this.ratios = ratios;
        this.oris = oris;
    }
    
    public void setArpSpanInfo(int[] mapQ, int[] weight, String[] itemTypes){
        arpSpanBpMapQ = mapQ;
        arpSpanBpWeight = weight;
        
        if (linkType == 1){
            arpSpanBpItem = itemTypes;
            if (itemTypes[0].contains("ARP") && !itemTypes[1].contains("ARP")){
                isPassed = false;
            }
            else if (!itemTypes[0].contains("ARP") && itemTypes[1].contains("APR")){
                isPassed = false;
            }
            else if (itemTypes[0].equals(itemTypes[1]) && !itemTypes[0].contains("ARP")){
                isPassed = false;
            }
        }            
        
    }
    
    public void setSelfLinkedInfo(int[] mapQ, int[] weight, String[] itemTypes){
        selfLinkedMapQ = mapQ;
        selfLinkedWeight = weight;
        
        if (linkType == -1){
            selfLinkedItemTypes = itemTypes;
            if (itemTypes[0].contains("ARP") && !itemTypes[1].contains("ARP")){
                isPassed = false;
            }
            else if (!itemTypes[0].contains("ARP") && itemTypes[1].contains("APR")){
                isPassed = false;
            }
            else if (itemTypes[0].equals(itemTypes[1]) && !itemTypes[0].contains("ARP")){
                isPassed = false;
            }
        }           
        
    }
    /**
     * Used to check if two SV are identical, just for SV < readLen
     * @param svInfo
     * @return 
     */
    public boolean identicalSV(svOutInfo svInfo){
        boolean identical = false;
        if (svInfo.start == this.start && svInfo.end == this.end){
            identical = true;
        }
        if (svInfo.end - svInfo.start == this.end - this.start){
            // sv shift to right by 1bp
            if (svInfo.start == this.start + 1 && svInfo.end == this.end + 1){
                identical = true;
            }
            // sv shift to left by 1bp
            if (svInfo.start == this.start - 1 && svInfo.end == this.end - 1){
                identical = true;
            }
        }
        return identical;
    }
       
    
  
    
    public void writeVariantsOutput(BufferedWriter regionWriter, String chrName, StringBuilder sb) throws IOException{
        if (end - start > 50){            
            if(linkType == 1 && (arpSpanBpMapQ[0] > arpSpanBpWeight[0] * 20 && arpSpanBpMapQ[1] > arpSpanBpWeight[1] * 20)){
                sb.append(chrName);
                sb.append(toString());
                regionWriter.write(sb.toString());
                regionWriter.newLine();
            }            
            else if (linkType != 1){
                sb.append(chrName);
                sb.append(toString());
                regionWriter.write(sb.toString());
                regionWriter.newLine();
            }            
        }
        if (end == -1){
            sb.append(chrName);
            sb.append(toString());
            regionWriter.write(sb.toString());
            regionWriter.newLine();
        }       
    }
}
