import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

public class SV_Candidate_Filter {

    public static void filter_sv_candidate(List list, String path) throws FileNotFoundException {

        List<SV_CandidateDeletion> candidateDeletions= (List <SV_CandidateDeletion>) list.get(0);
        List<SV_CandidateInsertion> candidateInsertions= (List <SV_CandidateInsertion>) list.get(1);
        List<SV_CandidateInversion> candidateInversions= (List <SV_CandidateInversion>) list.get(2);
        List<SV_CandidateDuplicationTandem> candidateDuplicationTandems= (List <SV_CandidateDuplicationTandem>) list.get(3);
        List<SV_CandidateInterspersedDuplication> candidateInterspersedDuplications= (List <SV_CandidateInterspersedDuplication>) list.get(4);
        List<SV_CandidateTranslocation> candidateTranslocations= (List <SV_CandidateTranslocation>) list.get(5);

        File deletion =new File(path+"_deletions.csv");
        File novelinsertion =new File(path+"_novelinsertions.csv");
        File inversion =new File(path+"_inversions.csv");
        File tandemduplication =new File(path+"_tandemduplications.csv");
        File int_duplication =new File(path+"_interspersedduplications.csv");
        File translocation =new File(path+"_translocations.csv");
        File block_interchange =new File(path+"_blockinterchanges.csv");

        PrintStream outputdeletion = new PrintStream(new FileOutputStream(deletion));
        PrintStream outputnovelinsertion = new PrintStream(new FileOutputStream(novelinsertion));
        PrintStream outputinversion = new PrintStream(new FileOutputStream(inversion));
        PrintStream outputtandemduplication = new PrintStream(new FileOutputStream(tandemduplication));
        PrintStream outputintduplication= new PrintStream(new FileOutputStream(int_duplication));
        PrintStream outputtranslocation= new PrintStream(new FileOutputStream(translocation));
        PrintStream outputblockinterchange= new PrintStream(new FileOutputStream(block_interchange));
        int candidateInversions_size=candidateInversions.size();
        outputinversion.println("contig"+"\t"+"start"+"\t"+"end"+"\t"+"score");
        for(int n=0;n<candidateInversions_size;n++){

            SV_CandidateInversion Inversion=candidateInversions.get(n);
            outputinversion.println(Inversion.contig+"\t"+String.valueOf(Inversion.start_position_on_ref)+"\t"+String.valueOf(Inversion.end_position_on_ref)+"\t"+String.valueOf(Inversion.score));

        }

        for(int i=0;i<candidateTranslocations.size();i++){

            SV_CandidateTranslocation trans=candidateTranslocations.get(i);
            String trans_source_contig=trans.contig1;
            long trans_source_start=trans.start_position_on_ref;
            long trans_source_end=trans.end_position_on_ref;
            String trans_destination_contig=trans.contig2;
            long trans_destination_start=trans.destination_start_position_on_ref;
            long trans_destination_end=trans.destination_end_position_on_ref;
            double trans_score=trans.score;

            for(int j=0;j<candidateTranslocations.size();j++){

                if(j!=i){

                    SV_CandidateTranslocation trans_1=candidateTranslocations.get(j);
                    String trans_source_contig_1=trans_1.contig1;
                    long trans_source_start_1=trans_1.start_position_on_ref;
                    long trans_source_end_1=trans_1.end_position_on_ref;
                    String trans_destination_contig_1=trans_1.contig2;
                    long trans_destination_start_1=trans_1.destination_start_position_on_ref;
                    long trans_destination_end_1=trans_1.destination_end_position_on_ref;
                    double trans_score_1=trans_1.score;

                    if(trans_source_contig.equals(trans_destination_contig)&&trans_source_contig_1.equals(trans_destination_contig_1)&&trans_source_contig.equals(trans_source_contig_1)){

                        if(Math.abs(trans_source_end-trans_destination_start_1)<=50&&Math.abs(trans_destination_start-trans_destination_end_1)<=50&&Math.abs(trans_destination_end-trans_source_start_1)<=50){

                            candidateInterspersedDuplications.add(new SV_CandidateInterspersedDuplication(trans_source_contig,trans_source_start,trans_source_end,trans_destination_contig,trans_destination_start,trans_destination_end,trans_score,false));
                            candidateInterspersedDuplications.add(new SV_CandidateInterspersedDuplication(trans_source_contig_1,trans_source_start_1,trans_source_end_1,trans_destination_contig_1,trans_destination_start_1,trans_destination_end_1,trans_score_1,false));

                            if(j>i){

                                candidateTranslocations.remove(j);
                                candidateTranslocations.remove(i);
                                i--;
                                break;

                            }
                            else{

                                candidateTranslocations.remove(i);
                                candidateTranslocations.remove(j);
                                i--;
                                break;

                            }

                        }

                    }

                }

            }

        }

        int duplication_tandem_signatures_size=candidateDuplicationTandems.size();

        for(int i=0;i<candidateInterspersedDuplications.size();i++){

            SV_CandidateInterspersedDuplication interspersedduplications=candidateInterspersedDuplications.get(i);
            String interspersedduplications_destination_contig=interspersedduplications.contig2;
            long interspersedduplications_destination_start=interspersedduplications.destination_start_position_on_ref;
            long interspersedduplications_destination_end=interspersedduplications.destination_end_position_on_ref;


            for(int j=0;j<duplication_tandem_signatures_size;j++){

                SV_CandidateDuplicationTandem duplicationTandem=candidateDuplicationTandems.get(j);
                String source_contig=duplicationTandem.contig;
                long source_start=duplicationTandem.start_position_on_ref;
                long source_end=duplicationTandem.end_position_on_ref;

                if(interspersedduplications_destination_contig.equals(source_contig)&&Math.abs(interspersedduplications_destination_start-source_end)<=50&&Math.abs(interspersedduplications_destination_end-source_start)<=50){

                    candidateInterspersedDuplications.remove(i);
                    i--;
                    break;

                }

            }

        }

        int interspersedduplicationsignatures_size=candidateInterspersedDuplications.size();

        for(int i=0;i<candidateDuplicationTandems.size();i++){

            SV_CandidateDuplicationTandem duplicationTandem=candidateDuplicationTandems.get(i);
            String source_contig=duplicationTandem.contig;
            long source_start=duplicationTandem.start_position_on_ref;
            long source_end=duplicationTandem.end_position_on_ref;
            boolean mark=false;

            for(int j=0;j<interspersedduplicationsignatures_size;j++){

                SV_CandidateInterspersedDuplication interspersedduplications=candidateInterspersedDuplications.get(j);
                String interspersedduplications_source_contig=interspersedduplications.contig1;
                long interspersedduplications_source_start=interspersedduplications.start_position_on_ref;
                long interspersedduplications_source_end=interspersedduplications.end_position_on_ref;
                long interspersedduplications_destination_start=interspersedduplications.destination_start_position_on_ref;
                long interspersedduplications_destination_end=interspersedduplications.destination_end_position_on_ref;

                if(source_contig.equals(interspersedduplications_source_contig)){

                    if(Math.abs(source_start-interspersedduplications_source_start)<=50&&Math.abs(source_end-interspersedduplications_destination_start)<=50){

                        candidateDuplicationTandems.remove(i);
                        i--;
                        mark=true;
                        break;

                    }
                    else if(Math.abs(source_start-interspersedduplications_destination_end)<=50&&Math.abs(source_end-interspersedduplications_source_end)<=50){

                        candidateDuplicationTandems.remove(i);
                        i--;
                        mark=true;
                        break;

                    }

                }

            }

            if(!mark){

                int size=candidateTranslocations.size();

                for(int j=0;j<size;j++){

                    SV_CandidateTranslocation translocation1=candidateTranslocations.get(j);
                    String translocation_source_contig=translocation1.contig1;
                    long translocation_source_start=translocation1.start_position_on_ref;
                    long translocation_source_end=translocation1.end_position_on_ref;
                    String translocation_destination_contig=translocation1.contig2;
                    long translocation_destination_start=translocation1.destination_start_position_on_ref;
                    long translocation_destination_end=translocation1.destination_end_position_on_ref;

                    if(source_contig.equals(translocation_source_contig)&&source_contig.equals(translocation_destination_contig)){

                        if(Math.abs(source_start-translocation_source_start)<=50&&Math.abs(source_end-translocation_destination_start)<=50){

                            candidateDuplicationTandems.remove(i);
                            i--;
                            break;

                        }
                        else if(Math.abs(source_start-translocation_destination_end)<=50&&Math.abs(source_end-translocation_source_end)<=50){

                            candidateDuplicationTandems.remove(i);
                            i--;
                            break;

                        }

                    }

                }

            }

        }

        SV_CandidateInterspersedDuplication current_int_duplication = null;
        SV_CandidateDuplicationTandem current_tan_duplication=null;
        duplication_tandem_signatures_size=candidateDuplicationTandems.size();

        for(int n=0;n<candidateInsertions.size();n++){

            SV_CandidateInsertion inserted_rigion=candidateInsertions.get(n);
            String contig1=inserted_rigion.contig;
            long start1=inserted_rigion.start_position_on_ref;
            long end1=inserted_rigion.end_position_on_ref;
            long length1=end1-start1;
            String contig2 = null;
            long start2=0;
            long end2 = 0;
            long length2;
            boolean mark=false;

            for(int i=0;i<interspersedduplicationsignatures_size;i++){

                current_int_duplication=candidateInterspersedDuplications.get(i);
                contig2=current_int_duplication.contig2;
                start2=current_int_duplication.destination_start_position_on_ref;
                end2=current_int_duplication.destination_end_position_on_ref;

                if(contig2.compareTo(contig1)<0||(contig2.equals(contig1)&&end2<start1)||(contig2.equals(contig1)&&end1<start2)){

                    continue;

                }

                length2=end2-start2;

                if(contig2.equals(contig1)&&start2<end1&&(length1-length2)/Math.max(length2,length1)<0.2){

                    candidateInsertions.remove(n);
                    n--;
                    mark=true;
                    break;

                }

            }

            if(!mark){

                for(int i=0;i<duplication_tandem_signatures_size;i++){

                    current_tan_duplication=candidateDuplicationTandems.get(i);
                    contig2=current_tan_duplication.contig;
                    start2=current_tan_duplication.end_position_on_ref;
                    end2=start2+(start2-current_tan_duplication.start_position_on_ref)*current_tan_duplication.copies;

                    if(contig2.compareTo(contig1)<0||(contig2.equals(contig1)&&end2<start1)||(contig2.equals(contig1)&&end1<start2)){

                        continue;

                    }

                    length2=end2-start2;

                    if(contig2.equals(contig1)&&start2<end1&&(length1-length2)/Math.max(length2,length1)<0.2){

                        candidateInsertions.remove(n);
                        n--;
                        break;

                    }

                }

            }

        }

        for(int n=0;n<candidateTranslocations.size();n++){

            SV_CandidateTranslocation check_translocation=candidateTranslocations.get(n);
            String translocation_source_contig=check_translocation.contig1;
            String translocation_destination_contig=check_translocation.contig2;
            long translocation_source_start=check_translocation.start_position_on_ref;
            long translocation_destination_start=check_translocation.destination_start_position_on_ref;
            long translocation_source_end=check_translocation.end_position_on_ref;
            long translocation_destination_end=check_translocation.destination_end_position_on_ref;
            String source_contig = null;
            long source_start=0;
            long source_end = 0 ;
            String destination_contig;
            long destination_start=0;
            long destination_end=0;
            boolean mark=false;

            for(int i=0;i<interspersedduplicationsignatures_size;i++){

                current_int_duplication=candidateInterspersedDuplications.get(i);
                source_contig=current_int_duplication.contig1;
                source_start=current_int_duplication.start_position_on_ref;
                source_end=current_int_duplication.end_position_on_ref;
                destination_contig=current_int_duplication.contig2;
                destination_start=current_int_duplication.destination_start_position_on_ref;
                destination_end=current_int_duplication.destination_end_position_on_ref;

                if( Math.abs(translocation_source_start-source_start)<=50&&Math.abs(translocation_source_end-source_end)<=50){

                    if(translocation_source_contig.equals(source_contig)){

                        candidateTranslocations.remove(n);
                        n--;
                        mark=true;
                        break;

                    }
                }
                else if((Math.abs(translocation_source_start-destination_end)<=50&&Math.abs(translocation_destination_start-source_end)<=50)||(Math.abs(translocation_source_end-destination_start)<=50&&Math.abs(translocation_destination_end-source_start)<=50)) {

                    if(translocation_destination_contig.equals(source_contig)&&translocation_source_contig.equals(destination_contig)){

                       candidateTranslocations.remove(n);
                       n--;
                       mark=true;
                       break;

                   }

                }

            }

            if(!mark){

                for(int i=0;i<duplication_tandem_signatures_size;i++){

                    current_tan_duplication=candidateDuplicationTandems.get(i);
                    source_contig=current_tan_duplication.contig;
                    source_start=current_tan_duplication.start_position_on_ref;
                    source_end=current_tan_duplication.end_position_on_ref;

                    if(translocation_source_contig.equals(source_contig)&&translocation_destination_contig.equals(source_contig)&&Math.abs(translocation_source_start-source_start)<=50&&Math.abs(translocation_source_end-source_end)<=50){

                        candidateTranslocations.remove(n);
                        n--;
                        break;

                    }

                }

            }

        }

        for(int i=0;i<candidateDeletions.size();i++){

            SV_CandidateDeletion check_deletion=candidateDeletions.get(i);
            String deletion_contig=check_deletion.contig;
            long deletion_start=check_deletion.start_position_on_ref;
            long deletion_end=check_deletion.end_position_on_ref;
            int size=candidateInterspersedDuplications.size();
            boolean mark=false;

            for(int j=0;j<size;j++){

                SV_CandidateInterspersedDuplication inters=candidateInterspersedDuplications.get(j);
                String source_contig=inters.contig1;
                long source_start=inters.start_position_on_ref;
                long source_end=inters.end_position_on_ref;
                long destination_start=inters.destination_start_position_on_ref;
                long destination_end=inters.destination_end_position_on_ref;

                if(deletion_contig.equals(source_contig)){

                    if(Math.abs(deletion_start-destination_start)<=50&&Math.abs(deletion_end-source_start)<=50){

                        candidateDeletions.remove(i);
                        i--;
                        mark=true;
                        break;

                    }
                    else if(Math.abs(deletion_start-source_end)<=50&&Math.abs(deletion_end-destination_end)<=50){

                        candidateDeletions.remove(i);
                        i--;
                        mark=true;
                        break;

                    }

                }

            }

            if(!mark){

                size=candidateTranslocations.size();

                for(int j=0;j<size;j++){

                    SV_CandidateTranslocation trans=candidateTranslocations.get(j);
                    String source_contig=trans.contig1;
                    long source_start=trans.start_position_on_ref;
                    long source_end=trans.end_position_on_ref;
                    String destination_contig=trans.contig2;
                    long destination_start=trans.destination_start_position_on_ref;
                    long destination_end=trans.destination_end_position_on_ref;

                    if(deletion_contig.equals(source_contig)&&deletion_contig.equals(destination_contig)){

                        if(Math.abs(deletion_start-destination_start)<=50&&Math.abs(deletion_end-source_start)<=50){

                            candidateDeletions.remove(i);
                            i--;
                            break;

                        }
                        else if(Math.abs(deletion_start-source_end)<=50&&Math.abs(deletion_end-destination_end)<=50){

                            candidateDeletions.remove(i);
                            i--;
                            break;

                        }

                    }

                }

            }

        }

        ArrayList<ArrayList> destinations_infor=new ArrayList <>();
        int int_dup_size=candidateInterspersedDuplications.size();

        for(int i=0;i<int_dup_size;i++){

            ArrayList destination_infor=new ArrayList();
            SV_CandidateInterspersedDuplication int_dup=candidateInterspersedDuplications.get(i);
            destination_infor.add(int_dup.contig2);
            destination_infor.add(int_dup.destination_start_position_on_ref);
            destination_infor.add(int_dup.destination_end_position_on_ref);
            destinations_infor.add(destination_infor);

        }

        int trans_size=candidateTranslocations.size();

        for(int i=0;i<trans_size;i++){

            ArrayList destination_infor=new ArrayList();
            SV_CandidateTranslocation trans=candidateTranslocations.get(i);
            destination_infor.add(trans.contig2);
            destination_infor.add(trans.destination_start_position_on_ref);
            destination_infor.add(trans.destination_end_position_on_ref);
            destinations_infor.add(destination_infor);

        }

        int candidateDeletions_size=candidateDeletions.size();
        boolean mark;
        for(int i=0;i<int_dup_size;i++){

            SV_CandidateInterspersedDuplication int_dup=candidateInterspersedDuplications.get(i);
            String int_source_contig=int_dup.contig1;
            long int_source_start=int_dup.start_position_on_ref;
            long int_source_end=int_dup.end_position_on_ref;
            mark=false;

            for(int n=0;n<candidateDeletions_size;n++){

                SV_CandidateDeletion del_candidate=candidateDeletions.get(n);
                String del_contig=del_candidate.contig;
                long del_start=del_candidate.start_position_on_ref;
                long del_end=del_candidate.end_position_on_ref;

                if(int_source_contig.equals(del_contig)&&Math.abs(int_source_start-del_start)<=50&&Math.abs(int_source_end-del_end)<=50){

                    candidateInterspersedDuplications.set(i,new SV_CandidateInterspersedDuplication(int_dup.contig1,int_dup.start_position_on_ref,int_dup.end_position_on_ref,int_dup.contig2,int_dup.destination_start_position_on_ref,int_dup.destination_end_position_on_ref,int_dup.score,true));
                    mark=true;
                    break;

                }

            }

            if(!mark){

                int destinations_infor_size=destinations_infor.size();

                for(int j=0;j<destinations_infor_size;j++){

                    ArrayList destination_infor=destinations_infor.get(j);

                    if(int_dup.contig1.equals(destination_infor.get(0))){

                        if((long)destination_infor.get(1)<=int_dup.start_position_on_ref&&int_dup.end_position_on_ref<=(long)destination_infor.get(2)){

                            candidateInterspersedDuplications.set(i,new SV_CandidateInterspersedDuplication(int_dup.contig1,int_dup.start_position_on_ref,int_dup.end_position_on_ref,int_dup.contig2,int_dup.destination_start_position_on_ref,int_dup.destination_end_position_on_ref,int_dup.score,true));
                            break;
//                               Global.interspersedduplicationsignatures.add(new SignatureInterspersedDuplication(Potential_Trans_Inters.contig1,Potential_Trans_Inters.start_position_on_ref,Potential_Trans_Inters.end_position_on_ref,Potential_Trans_Inters.contig2,Potential_Trans_Inters.destination_start_position_on_ref,Potential_Trans_Inters.destination_end_position_on_ref,Potential_Trans_Inters.signature,"cut&paste",Potential_Trans_Inters.confidence_level));

                        }

                    }

                }

            }

        }

        for(int i=0;i<candidateInterspersedDuplications.size()-1;i++){

            SV_CandidateInterspersedDuplication is_cut_paste=candidateInterspersedDuplications.get(i);

            if(is_cut_paste.is_cutpaste){

                String is_cut_paste_source_contig=is_cut_paste.contig1;
                long is_cut_paste_source_start=is_cut_paste.start_position_on_ref;
                long is_cut_paste_source_end=is_cut_paste.end_position_on_ref;
                String is_cut_paste_destination_contig=is_cut_paste.contig2;
                long is_cut_paste_destination_start=is_cut_paste.destination_start_position_on_ref;
                long is_cut_paste_destination_end=is_cut_paste.destination_end_position_on_ref;
                double  is_cut_paste_score=is_cut_paste.score;

                for(int j=i+1;j<candidateInterspersedDuplications.size();j++){

                    SV_CandidateInterspersedDuplication is_cut_paste_1=candidateInterspersedDuplications.get(j);

                    if(is_cut_paste_1.is_cutpaste){

                        String is_cut_paste_1_source_contig=is_cut_paste_1.contig1;
                        long is_cut_paste_1_source_start=is_cut_paste_1.start_position_on_ref;
                        long is_cut_paste_1_source_end=is_cut_paste_1.end_position_on_ref;
                        String is_cut_paste_1_destination_contig=is_cut_paste_1.contig2;
                        long is_cut_paste_1_destination_start=is_cut_paste_1.destination_start_position_on_ref;
                        long is_cut_paste_1_destination_end=is_cut_paste_1.destination_end_position_on_ref;
                        double  is_cut_paste_1_score=is_cut_paste_1.score;

                        if(is_cut_paste_source_contig.equals(is_cut_paste_1_destination_contig)&&is_cut_paste_destination_contig.equals(is_cut_paste_1_source_contig)){

                            if(Math.abs(is_cut_paste_destination_start-is_cut_paste_1_source_start)<=50&&Math.abs(is_cut_paste_1_destination_start-is_cut_paste_source_start)<=50){

                                double score=(is_cut_paste_score+is_cut_paste_1_score)/2;
                                candidateTranslocations.add(new SV_CandidateTranslocation(is_cut_paste_source_contig,is_cut_paste_source_start,is_cut_paste_source_end,is_cut_paste_1_source_contig,is_cut_paste_1_source_start,is_cut_paste_1_source_end,score));
                                candidateInterspersedDuplications.remove(j);
                                candidateInterspersedDuplications.remove(i);
                                i--;
                                break;

                            }

                        }

                    }

                }

            }

        }

        outputdeletion.println("contig"+"\t"+"start"+"\t"+"end"+"\t"+"score"+"\t"+"span");

        for(int n=0;n<candidateDeletions_size;n++){

            SV_CandidateDeletion Deletion=candidateDeletions.get(n);
            outputdeletion.println(Deletion.contig+"\t"+String.valueOf(Deletion.start_position_on_ref)+"\t"+String.valueOf(Deletion.end_position_on_ref)+"\t"+String.valueOf(Deletion.score)+"\t"+String.valueOf(Deletion.end_position_on_ref-Deletion.start_position_on_ref));

        }

        int candidateInsertions_size=candidateInsertions.size();
        outputnovelinsertion.println("contig"+"\t"+"start"+"\t"+"end"+"\t"+"score"+"\t"+"span");

        for(int n=0;n<candidateInsertions_size;n++){

            SV_CandidateInsertion Insertion=candidateInsertions.get(n);
            outputnovelinsertion.println(Insertion.contig+"\t"+String.valueOf(Insertion.start_position_on_ref)+"\t"+String.valueOf(Insertion.end_position_on_ref)+"\t"+String.valueOf(Insertion.score)+"\t"+String.valueOf(Insertion.end_position_on_ref-Insertion.start_position_on_ref));

        }

        duplication_tandem_signatures_size=candidateDuplicationTandems.size();
        outputtandemduplication.println("contig"+"\t"+"start"+"\t"+"end"+"\t"+"copies"+"\t"+"score");

        for(int n=0;n<duplication_tandem_signatures_size;n++){

            SV_CandidateDuplicationTandem tandem_duplication=candidateDuplicationTandems.get(n);
            outputtandemduplication.println(tandem_duplication.contig+"\t"+String.valueOf(tandem_duplication.start_position_on_ref)+"\t"+String.valueOf(tandem_duplication.end_position_on_ref)+"\t"+String.valueOf(tandem_duplication.copies)+"\t"+String.valueOf(tandem_duplication.score));

        }

        interspersedduplicationsignatures_size=candidateInterspersedDuplications.size();
        outputintduplication.println("contig"+"\t"+"source_start"+"\t"+"source_end"+"\t"+"dest_contig"+"\t"+"dest_start"+"\t"+"dest_end"+"\t"+"score"+"\t"+"cut&paste");
        for(int n=0;n<interspersedduplicationsignatures_size;n++){

            SV_CandidateInterspersedDuplication DuplicationInterspersed=candidateInterspersedDuplications.get(n);
            outputintduplication.println(DuplicationInterspersed.contig1+"\t"+String.valueOf(DuplicationInterspersed.start_position_on_ref)+"\t"+String.valueOf(DuplicationInterspersed.end_position_on_ref)+"\t"+DuplicationInterspersed.contig2+"\t"+String.valueOf(DuplicationInterspersed.destination_start_position_on_ref)+"\t"+String.valueOf(DuplicationInterspersed.destination_end_position_on_ref)+"\t"+String.valueOf(DuplicationInterspersed.score)+"\t"+String.valueOf(DuplicationInterspersed.is_cutpaste));

        }

        int candidateTranslocations_size=candidateTranslocations.size();
        outputtranslocation.println("contig"+"\t"+"source_start"+"\t"+"source_end"+"\t"+"dest_contig"+"\t"+"dest_start"+"\t"+"dest_end"+"\t"+"score");
        outputblockinterchange.println("contig"+"\t"+"source_start"+"\t"+"source_end"+"\t"+"dest_contig"+"\t"+"dest_start"+"\t"+"dest_end"+"\t"+"score");

        for(int n=0;n<candidateTranslocations_size;n++){

            SV_CandidateTranslocation Translocation=candidateTranslocations.get(n);

            if(Translocation.contig1.equals(Translocation.contig2)){

                outputblockinterchange.println(Translocation.contig1+"\t"+String.valueOf(Translocation.start_position_on_ref)+"\t"+String.valueOf(Translocation.end_position_on_ref)+"\t"+Translocation.contig2+"\t"+String.valueOf(Translocation.destination_start_position_on_ref)+"\t"+String.valueOf(Translocation.destination_end_position_on_ref)+"\t"+String.valueOf(Translocation.score));

            }
            else{

                outputtranslocation.println(Translocation.contig1+"\t"+String.valueOf(Translocation.start_position_on_ref)+"\t"+String.valueOf(Translocation.end_position_on_ref)+"\t"+Translocation.contig2+"\t"+String.valueOf(Translocation.destination_start_position_on_ref)+"\t"+String.valueOf(Translocation.destination_end_position_on_ref)+"\t"+String.valueOf(Translocation.score));

            }

        }

    }

    private static class sort_distance implements Comparator<ArrayList> {

        public int compare(ArrayList o1, ArrayList o2) {

            if((double)o1.get(1)>(double)o2.get(1)){

                return 1;
            }
            else if((double)o1.get(1)==(double)o2.get(1)){


                return 0;
            }
            else{

                return -1;

            }

        }

    }

}
