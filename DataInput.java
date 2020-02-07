import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import java.io.*;
import java.sql.Time;
import java.util.*;

class parse_mapped_data implements Serializable {

    private long Score;
    private double mismap;
    private String Ref;
    private String Read;
    private long RefStartC;
    private long ReadStartC;
    private long RefAlignLen;
    private long ReadAlignLen;
    private long RefEndC;
    private long ReadEndC;
    private String RefStrand;
    private String ReadStrand;
    private long ReadLen;
    private String RefSeq;
    private String ReadSeq;
    private double Identity;
    private String CIGAR;
//    private double confidence_level=0;
//
//    public void setConfidence_level(double confidence_level) {
//
//        this.confidence_level = confidence_level;
//
//    }
//
//    public double getConfidence_level() {
//
//        return confidence_level;
//
//    }

    public void setCIGAR(String CIGAR) {

        this.CIGAR = CIGAR;
    }

    public String getCIGAR() {

        return CIGAR;

    }

    public void setScore(long Score) {
        this.Score = Score;
    }

    public long getScore() {
        return Score;
    }

    public void setMismap(double mismap) {
        this.mismap = mismap;
    }

    public double getMismap() {
        return mismap;
    }

    public void setRef(String Ref) {
        this.Ref = Ref;
    }

    public String getRef() {
        return Ref;
    }

    public void setRead(String Read) {
        this.Read = Read;
    }

    public String getRead() {
        return Read;
    }

    public void setRefStartC(long RefStartC) {
        this.RefStartC = RefStartC;
    }

    public long getRefStartC() {
        return RefStartC;
    }

    public void setReadStartC(long ReadStartC) {
        this.ReadStartC = ReadStartC;
    }

    public long getReadStartC() {
        return ReadStartC;
    }

    public void setRefAlignLen(long RefAlignLen) {
        this.RefAlignLen = RefAlignLen;
    }

    public long getRefAlignLen() {
        return RefAlignLen;
    }

    public void setReadAlignLen(long ReadAlignLen) {
        this.ReadAlignLen = ReadAlignLen;
    }

    public long getReadAlignLen() {
        return ReadAlignLen;
    }

    public void setRefEndC(long RefEndC) {
        this.RefEndC = RefEndC;
    }

    public long getRefEndC() {
        return RefEndC;
    }

    public void setReadEndC(long ReadEndC) {
        this.ReadEndC = ReadEndC;
    }

    public long getReadEndC() {
        return ReadEndC;
    }

    public void setRefStrand(String RefStrand) {
        this.RefStrand = RefStrand;
    }

    public String getRefStrand() {
        return RefStrand;
    }

    public void setReadStrand(String ReadStrand) {
        this.ReadStrand = ReadStrand;
    }

    public String getReadStrand() {
        return ReadStrand;
    }

    public void setReadLen(long Readlen) {
        this.ReadLen = Readlen;
    }

    public long getReadLen() {
        return ReadLen;
    }

    public void setRefSeq(String RefSeq) {
        this.RefSeq = RefSeq;
    }

    public String getRefSeq() {
        return RefSeq;
    }

    public void setReadSeq(String ReadSeq) {
        this.ReadSeq = ReadSeq;
    }

    public String getReadSeq() {
        return ReadSeq;
    }

    public void setIdentity(double Identity) {
        this.Identity = Identity;
    }

    public double getIdentity() {
        return Identity;
    }
}

class CorOrder implements Comparator<parse_mapped_data> {

    @Override
    public int compare(parse_mapped_data o1, parse_mapped_data o2) {

        if(o1.getReadStartC() > o2.getReadStartC()){

            return 1;

        }
        else if(o1.getReadStartC() == o2.getReadStartC()&&o1.getReadEndC() > o2.getReadEndC()){

            return 1;

        }
        else if(o1.getReadStartC() == o2.getReadStartC()&&o1.getReadEndC() == o2.getReadEndC()){

            return 0;

        }
        else {

            return -1;

        }

    }

}

class FilterThreshold{

    private static double MinIdentity=55;
    public double getMinIdentity() {
        return MinIdentity ;
    }

}


public class DataInput {

    public static void main(String [] args) throws IOException, ClassNotFoundException {

        Global.deletion_signatures=new LinkedList<>();
        Global.insertion_signatures=new LinkedList <>();
        Global.inversion_signatures=new LinkedList <>();
        Global.duplication_tandem_signatures=new LinkedList <>();
        Global.TranslocationBreakpoint_signatures=new LinkedList <>();
        Global.Potential_Trans_Inters_signatures=new LinkedList <>();
        Parameter_Setting Parameter_Setting=new Parameter_Setting();
//        Scanner datainput=new Scanner(System.in);
        long startTime = System.currentTimeMillis();
        String path=args[0];
        String path_split[]=path.split("\\.");
        String path_format=path_split[path_split.length-1];
        if(path_format.equals("maf")){

            alignDataInput_maf(path,Parameter_Setting);

        }
        else if(path_format.equals("sam")){

            alignDataInput_sam(path,Parameter_Setting);

        }
        SV_Candidate_Filter.filter_sv_candidate(T_I_C_SVClusters_Distinguish.parse_Potential_Trans_Inters_clusters(SV_Signature_Cluster.cluster_sv_signature(Parameter_Setting),Parameter_Setting),path_split[0]);
        long endTime = System.currentTimeMillis();
        System.out.println("程序运行时间：" + (endTime - startTime) + "ms");

    }

    public static void alignDataInput_maf(String path, Parameter_Setting Parameter_Setting) throws IOException {

        File MAFfile = new File(path);
        BufferedReader MAFbufferedReader;
        MAFbufferedReader = new BufferedReader(new FileReader(MAFfile));
        FilterThreshold filterthreshold = new FilterThreshold();
        String          MAFAlignline;
        String          Aline;
        String          Sline1;
        String          Sline2;
        String          currentRead = null;
        List<parse_mapped_data> currentReadStore=new LinkedList <>();
        String          [] charTest;
        String          [] midStore;
        long            refendc;
        long            readendc;
        long            matchcount;
        long            mismatchcount;
        long            insertioncount;
        long            deletioncount;
        long            alignlen;
        double          identity;
        parse_mapped_data AlignData[]=new parse_mapped_data[500];
        for(int i=0;i<500;i++){

            AlignData[i]=new parse_mapped_data();

        }
        long base_best_score_sum=0;
//        int j=0;
        int count = 1;
        int index=0;

        while((MAFAlignline = MAFbufferedReader.readLine()) != null){

            if(count <= 31){

                if(count==10){

                    String [] get_per_base_best_score=MAFAlignline.split("\\s+");
                    base_best_score_sum+=Long.parseLong(get_per_base_best_score[2]);

                }
                if(count==11){

                    String [] get_per_base_best_score=MAFAlignline.split("\\s+");
                    base_best_score_sum+=Long.parseLong(get_per_base_best_score[3]);

                }
                if(count==12){

                    String [] get_per_base_best_score=MAFAlignline.split("\\s+");
                    base_best_score_sum+=Long.parseLong(get_per_base_best_score[4]);

                }
                if(count==13){

                    String [] get_per_base_best_score=MAFAlignline.split("\\s+");
                    base_best_score_sum+=Long.parseLong(get_per_base_best_score[5]);

                }

                count++;
                continue;
            }

            Global.per_base_best_score=base_best_score_sum/4;
            matchcount = 0;
            mismatchcount = 0;
            insertioncount = 0;
            deletioncount = 0;
            charTest = MAFAlignline.split("\\s+");

            if(MAFAlignline.length() != 0&&!(charTest[0].equals("#"))){

                Aline = MAFAlignline.trim().replaceAll("="," ");
                midStore = Aline.split("\\s+");
                AlignData[index].setScore(Long.parseLong(midStore[2]));
                AlignData[index].setMismap(Double.parseDouble(midStore[4]));
                Sline1  =   MAFbufferedReader.readLine();
                midStore=   Sline1.split("\\s+");
                AlignData[index].setRef(midStore[1]);
                AlignData[index].setRefStartC(Long.parseLong(midStore[2]));
                AlignData[index].setRefAlignLen(Long.parseLong(midStore[3]));
                refendc =   AlignData[index].getRefStartC()+AlignData[index].getRefAlignLen()-1;
                AlignData[index].setRefEndC(refendc);
                AlignData[index].setRefStrand(midStore[4]);
                AlignData[index].setRefSeq(midStore[6].toUpperCase());
                Sline2  =   MAFbufferedReader.readLine();
                midStore=   Sline2.split("\\s+");
                AlignData[index].setRead(midStore[1]);
                AlignData[index].setReadStartC(Long.parseLong(midStore[2]));
                AlignData[index].setReadAlignLen(Long.parseLong(midStore[3]));
                readendc=   AlignData[index].getReadStartC()+AlignData[index].getReadAlignLen()-1;
                AlignData[index].setReadEndC(readendc);
                AlignData[index].setReadStrand(midStore[4]);
                AlignData[index].setReadLen(Long.parseLong(midStore[5]));
                AlignData[index].setReadSeq(midStore[6].toUpperCase());
                MAFbufferedReader.readLine();
                MAFbufferedReader.readLine();

                String CIGAR = null;
                String type=" ";
                String currtype=" ";
                long curtypecount=0;

                for (int i=0;i<AlignData[index].getRefSeq().length()-1&&i<AlignData[index].getReadSeq().length()-1;i++){

                    String Refbase=AlignData[index].getRefSeq().substring(i,i+1);
                    String querybase=AlignData[index].getReadSeq().substring(i,i+1);

                    if(Refbase.equals("-")){

                        if(querybase.equals("-")){


                        }
                        else{

                            currtype="I";
                            insertioncount++;

                        }

                    }
                    else{

                        if (querybase.equals("-")){

                            currtype="D";
                            deletioncount++;

                        }
                        else{

                            if(Refbase.equals(querybase)){

                                currtype="M";
                                matchcount++;

                            }else{

                                currtype="M";
                                mismatchcount++;

                            }

                        }

                    }

                    if(currtype.equals(type)){

                        curtypecount++;

                    }
                    else{

                        if(CIGAR==null){

                            CIGAR=curtypecount+type;
                        }
                        else{

                            CIGAR=CIGAR+curtypecount+type;

                        }
                        type=currtype;
                        curtypecount=1;

                    }

                }

                alignlen    =  matchcount+mismatchcount+insertioncount+deletioncount;
                identity    = (matchcount * 100.0) / alignlen;
                AlignData[index].setIdentity(identity);
                CIGAR=CIGAR+curtypecount+type;
                String  CIGARStore[]=CIGAR.split(" ");
                AlignData[index].setCIGAR(CIGARStore[1].trim());

                if(AlignData[index].getReadStrand().equals("-")){

                    AlignData[index].setRefStrand("-");
                    AlignData[index].setReadStrand("+");
                    AlignData[index].setReadStartC(AlignData[index].getReadLen()-1-AlignData[index].getReadEndC());
                    AlignData[index].setReadEndC(AlignData[index].getReadStartC()+AlignData[index].getReadAlignLen()-1);
//                    AlignData.setReadSeq(new StringBuffer(AlignData.getReadSeq()).reverse().toString().replaceAll("A","B").replaceAll("T","A").replaceAll("B","T").replaceAll("C","B").replaceAll("G","C").replaceAll("B","G"));

                }

                AlignData[index].setReadSeq(null);
                AlignData[index].setRefSeq(null);

                if(AlignData[index].getIdentity()>=filterthreshold.getMinIdentity()){

                    if(currentReadStore.size()!=0){

                        if(currentRead!=null&&!currentRead.equals(AlignData[index].getRead())){

                            if(currentReadStore.size()>7){

                                Collections.sort(currentReadStore,new Mismap());
                                currentReadStore=currentReadStore.subList(0,7);

                            }

                            Collections.sort(currentReadStore, new CorOrder());
//                            getCandidatemetrics result=getCandidatemetrics.getCandidatemetrics(currentReadStore);
//                            int currentReadStore_size=currentReadStore.size();
//                            double per_base_score=result.getPer_base_score();
//
//                            for(int m=0;m<currentReadStore_size;m++){
//
//                                currentReadStore.get(m).setConfidence_level(per_base_score);
//
//                            }

                            SV_Signature_Distinguish.analyze_read_segments(currentReadStore, Parameter_Setting);
//                            System.out.println("第"+String.valueOf(j++)+"条read处理完毕");
                            currentRead=AlignData[index].getRead();
                            currentReadStore=new LinkedList <>();
//                            AlignData[0]=AlignData[index];
//                            AlignData[0] = new parse_mapped_data(AlignData[index]);
                            copy(AlignData[0],AlignData[index]);
                            index=0;

                        }

                        currentReadStore.add(AlignData[index]);
                        index++;

                    }
                    else{

                        currentRead=AlignData[index].getRead();
                        currentReadStore.add(AlignData[index]);
                        index++;

                    }

                }

            }

        }

        if(MAFbufferedReader!=null){

            MAFbufferedReader.close();

        }

        if(currentReadStore!=null){

            if(currentReadStore.size()>7){

                Collections.sort(currentReadStore,new Mismap());
                currentReadStore=currentReadStore.subList(0,7);

            }

            Collections.sort(currentReadStore, new CorOrder());
//            getCandidatemetrics result=getCandidatemetrics.getCandidatemetrics(currentReadStore);
//            int currentReadStore_size=currentReadStore.size();
//            double per_base_score=result.getPer_base_score();
//
//            for(int m=0;m<currentReadStore_size;m++){
//
//                currentReadStore.get(m).setConfidence_level(per_base_score);
//
//            }

            SV_Signature_Distinguish.analyze_read_segments(currentReadStore, Parameter_Setting);
//            System.out.println("第"+String.valueOf(j++)+"条read处理完毕");

        }

    }

    private static void copy(parse_mapped_data copy, parse_mapped_data source) {

        copy.setScore(source.getScore());
        copy.setMismap(source.getMismap());
        copy.setRef(source.getRef());
        copy.setRefStartC(source.getRefStartC());
        copy.setRefAlignLen(source.getRefAlignLen());
        copy.setRefEndC(source.getRefEndC());
        copy.setRefStrand(source.getRefStrand());
        copy.setRefSeq(source.getRefSeq());
        copy.setRead(source.getRead());
        copy.setReadStartC(source.getReadStartC());
        copy.setReadAlignLen(source.getReadAlignLen());
        copy.setReadEndC(source.getReadEndC());
        copy.setReadStrand(source.getReadStrand());
        copy.setReadLen(source.getReadLen());
        copy.setReadSeq(source.getReadSeq());

    }

    public static void alignDataInput_sam(String path, Parameter_Setting Parameter_Setting) throws IOException {

        File SAMfile = new File(path);
        BufferedReader SAMbufferedReader;
        SAMbufferedReader = new BufferedReader(new FileReader(SAMfile));
        String          SAMAlignline;
        String          currentRead = null;
        List<parse_mapped_data> currentReadStore=new LinkedList <>();
        String          [] midStore;
        long            refendc;
        long            readendc;
        parse_mapped_data AlignData;
        int j=0;

        while((SAMAlignline = SAMbufferedReader.readLine()) != null){

            Global.per_base_best_score=2;
            AlignData = new parse_mapped_data();

            if(SAMAlignline.length() != 0&&!(SAMAlignline.substring(0,1).equals("@"))){

                midStore = SAMAlignline.trim().split("\\s+");

                if(midStore[1].equals("4")){

                    continue;

                }

                String Aliscore[]=midStore[11].split(":");
                AlignData.setScore(Long.parseLong(Aliscore[2]));
                AlignData.setRef(midStore[2]);
                AlignData.setRefStartC(Long.parseLong(midStore[3]));
                AlignData.setMismap(1/Double.parseDouble(midStore[4]));
                long RefAlignLen = 0;
                long ReadAlignLen = 0;
                String cigar=midStore[5];
                long length=cigar.length();
                int startpos=0;
                for(int i=0;i<length;i++){

                    if(i==0){

                        startpos=i;

                    }
                    else if('0'<=cigar.charAt(i)&&cigar.charAt(i)<='9'){

                        continue;

                    }
                    else if(cigar.charAt(i)=='S'&&i!=length-1){

                        AlignData.setReadStartC(Long.parseLong(cigar.substring(startpos,i)));
                        startpos=i+1;

                    }
                    else if(cigar.charAt(i)=='M'){

                        ReadAlignLen+=Long.parseLong(cigar.substring(startpos,i));
                        RefAlignLen+=Long.parseLong(cigar.substring(startpos,i));
                        startpos=i+1;

                    }
                    else if(cigar.charAt(i)=='D'){

                        RefAlignLen+=Long.parseLong(cigar.substring(startpos,i));
                        startpos=i+1;

                    }
                    else if(cigar.charAt(i)=='I'){

                        ReadAlignLen+=Long.parseLong(cigar.substring(startpos,i));
                        startpos=i+1;

                    }

                }

                AlignData.setRefAlignLen(RefAlignLen);
                refendc =   AlignData.getRefStartC()+AlignData.getRefAlignLen()-1;
                AlignData.setRefEndC(refendc);
                AlignData.setRefStrand("+");
                AlignData.setRead(midStore[0]);
                AlignData.setReadAlignLen(ReadAlignLen);
                readendc=   AlignData.getReadStartC()+AlignData.getReadAlignLen()-1;
                AlignData.setReadEndC(readendc);
                if(midStore[1].equals("0")){

                    AlignData.setReadStrand("+");

                }
                else if(midStore[1].equals("16")){

                    AlignData.setReadStrand("-");

                }
                else if(midStore[1].equals("2048")){

                    AlignData.setReadStrand("+");

                }
                else if(midStore[1].equals("2064")){

                    AlignData.setReadStrand("-");

                }
                AlignData.setReadLen(midStore[9].length());
                AlignData.setCIGAR(cigar);

                if(AlignData.getReadStrand().equals("-")){

                    AlignData.setRefStrand("-");
                    AlignData.setReadStrand("+");
                    AlignData.setReadStartC(AlignData.getReadLen()-1-AlignData.getReadEndC());
                    AlignData.setReadEndC(AlignData.getReadStartC()+AlignData.getReadAlignLen()-1);
//                    AlignData.setReadSeq(new StringBuffer(AlignData.getReadSeq()).reverse().toString().replaceAll("A","B").replaceAll("T","A").replaceAll("B","T").replaceAll("C","B").replaceAll("G","C").replaceAll("B","G"));

                }

                if(currentReadStore.size()!=0){

                    if(currentRead!=null&&!currentRead.equals(AlignData.getRead())){

                        if(currentReadStore.size()>7){

                            Collections.sort(currentReadStore,new Mismap());
                            currentReadStore=currentReadStore.subList(0,7);

                        }

                        Collections.sort(currentReadStore, new CorOrder());
//                        getCandidatemetrics result=getCandidatemetrics.getCandidatemetrics(currentReadStore);
//                        int currentReadStore_size=currentReadStore.size();
//                        double per_base_score=result.getPer_base_score();
//
//                        for(int m=0;m<currentReadStore_size;m++){
//
//                            currentReadStore.get(m).setConfidence_level(per_base_score);
//
//                        }

                        SV_Signature_Distinguish.analyze_read_segments(currentReadStore, Parameter_Setting);
                        System.out.println("第"+String.valueOf(j++)+"条read处理完毕");
                        currentRead=AlignData.getRead();
                        currentReadStore=new LinkedList <>();

                    }

                    currentReadStore.add(AlignData);

                }
                else{

                    currentRead=AlignData.getRead();
                    currentReadStore.add(AlignData);

                }

            }

        }

        if(SAMbufferedReader!=null){

            SAMbufferedReader.close();

        }

        if(currentReadStore!=null){

            if(currentReadStore.size()>7){

                Collections.sort(currentReadStore,new Mismap());
                currentReadStore=currentReadStore.subList(0,7);

            }

            Collections.sort(currentReadStore, new CorOrder());
//            getCandidatemetrics result=getCandidatemetrics.getCandidatemetrics(currentReadStore);
//            int currentReadStore_size=currentReadStore.size();
//            double per_base_score=result.getPer_base_score();
//
//            for(int m=0;m<currentReadStore_size;m++){
//
//                currentReadStore.get(m).setConfidence_level(per_base_score);
//
//            }

            SV_Signature_Distinguish.analyze_read_segments(currentReadStore, Parameter_Setting);
            System.out.println("第"+String.valueOf(j++)+"条read处理完毕");

        }

    }
    private static class Mismap implements Comparator<parse_mapped_data>{

        public int compare(parse_mapped_data o1, parse_mapped_data o2){

            if(o1.getMismap()>o2.getMismap()){

                return 1;

            }
            else  if(o1.getMismap()==o2.getMismap()){

                return 1;

            }
            else{

                return -1;

            }

        }

    }

}
