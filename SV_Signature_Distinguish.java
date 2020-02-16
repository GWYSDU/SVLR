import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

class potentialshiftsbreakpointreadCOrder implements Comparator<ArrayList> {

    @Override
    public int compare(ArrayList o1, ArrayList o2) {

        if((long)o1.get(5)>(long)o2.get(5)){

            return 1;

        }
        else if((long)o1.get(5)==(long)o2.get(5)){

            return 0;

        }
        else {

            return -1;

        }

    }

}

class potentialtranslocationsbreakpointreadCorOrder implements Comparator<ArrayList> {

    @Override
    public int compare(ArrayList o1, ArrayList o2) {

        if((long)o1.get(6)>(long)o2.get(6)){

            return 1;

        }
        else if((long)o1.get(6)==(long)o2.get(6)){

            return 0;

        }
        else {

            return -1;

        }

    }

}
public class SV_Signature_Distinguish{

    public static void analyze_read_segments(List <parse_mapped_data> parse_mapped_data, Parameter_Setting parameter_setting) throws FileNotFoundException {

        List<ArrayList> tandem_duplications=new LinkedList<>();
        List<ArrayList> potentialinsertions=new LinkedList <>();
        List<ArrayList> potentialdeletions=new LinkedList <>();
        List<ArrayList>potentialtranslocationsbreakpoint=new LinkedList <>();
        boolean is_used=false;
        if(parse_mapped_data!=null){

            if(parse_mapped_data.size()>0){

                String read_name=parse_mapped_data.get(0).getRead();

                for(int n=0;n<parse_mapped_data.size();n++){

                    parse_mapped_data parseMappedData=parse_mapped_data.get(n);
                    String cigar=parseMappedData.getCIGAR();
                    long length=cigar.length();
                    int startpos=0;
                    long RefAlignLen = 0;
                    for(int i=0;i<length;i++){

                        if(i==0){

                            startpos=i;

                        }
                        else if('0'<=cigar.charAt(i)&&cigar.charAt(i)<='9'){

                            continue;

                        }
                        else if(cigar.charAt(i)=='D'){

                            long gaplength=Long.parseLong(cigar.substring(startpos,i));

                            if(gaplength>=parameter_setting.min_sv_size){

                                ArrayList potentialdeletion=new ArrayList();
                                potentialdeletion.add(parseMappedData.getRef());
                                potentialdeletion.add(parseMappedData.getRefStartC()+RefAlignLen);
                                potentialdeletion.add(parseMappedData.getRefStartC()+RefAlignLen+gaplength);
                                potentialdeletion.add(gaplength);
                                potentialdeletion.add("cigar");
                                potentialdeletion.add("fwd");
                                potentialdeletion.add(((double)parseMappedData.getScore())/parseMappedData.getReadAlignLen());
                                potentialdeletions.add(potentialdeletion);

                            }

                            RefAlignLen+=gaplength;
                            startpos=i+1;

                        }
                        else if(cigar.charAt(i)=='I'){

                            long gaplength=Long.parseLong(cigar.substring(startpos,i));

                            if(gaplength>=parameter_setting.min_sv_size){

                                ArrayList potentialinsertion=new ArrayList();
                                potentialinsertion.add(parseMappedData.getRef());
                                potentialinsertion.add(parseMappedData.getRefStartC()+RefAlignLen);
                                potentialinsertion.add(parseMappedData.getRefStartC()+RefAlignLen);
                                potentialinsertion.add(Long.parseLong(cigar.substring(startpos,i)));
                                potentialinsertion.add("cigar");
                                potentialinsertion.add("fwd");
                                potentialinsertion.add(((double)parseMappedData.getScore())/parseMappedData.getReadAlignLen());
                                potentialinsertions.add(potentialinsertion);

                            }

                            startpos=i+1;

                        }
                        else if(cigar.charAt(i)!='S'&&cigar.charAt(i)!='H'){

                            RefAlignLen+=Long.parseLong(cigar.substring(startpos,i));
                            startpos=i+1;

                        }
                        else{

                            startpos=i+1;

                        }

                    }

                }

                for(int n=0;n<parse_mapped_data.size()-1;n++){

                    parse_mapped_data alignment_current=parse_mapped_data.get(n);
                    parse_mapped_data alignment_next=parse_mapped_data.get(n+1);
                    double alignment_current_confidence_level=((double)alignment_current.getScore())/alignment_current.getReadAlignLen();
                    double alignment_next_confidence_level=((double)alignment_next.getScore())/alignment_next.getReadAlignLen();
                    double confidence_level=(alignment_current_confidence_level+alignment_next_confidence_level)/2;
                    long distance_on_read=alignment_next.getReadStartC()-alignment_current.getReadEndC();
                    long distance_on_reference;
                    String ref_chr;

                    if(alignment_current.getRef().equals(alignment_next.getRef())){

                        ref_chr=alignment_current.getRef();

                        if(alignment_current.getRefStrand().equals(alignment_next.getRefStrand())){

                            if(alignment_current.getRefStrand().equals("-")){

                                distance_on_reference=alignment_current.getRefStartC()-alignment_next.getRefEndC();

                            }
                            else{

                                distance_on_reference=alignment_next.getRefStartC()-alignment_current.getRefEndC();

                            }

                            if(distance_on_read>= -parameter_setting.segment_overlap_tolerance){

                                if(distance_on_reference>= -parameter_setting.segment_overlap_tolerance){

                                    long deviation=distance_on_read-distance_on_reference;

                                    if(deviation>=parameter_setting.min_sv_size){

                                        if(distance_on_reference<=parameter_setting.segment_gap_tolerance){

                                            if(alignment_current.getRefStrand().equals("+")){

                                                ArrayList potentialinsertion=new ArrayList();
                                                potentialinsertion.add(ref_chr);
                                                potentialinsertion.add(alignment_current.getRefEndC());
                                                potentialinsertion.add(alignment_next.getRefStartC());
                                                potentialinsertion.add(Math.abs(deviation));
                                                potentialinsertion.add("suppl");
                                                potentialinsertion.add("fwd");
                                                potentialinsertion.add(confidence_level);
                                                potentialinsertions.add(potentialinsertion);

                                            }
                                            else{

                                                ArrayList potentialinsertion=new ArrayList();
                                                potentialinsertion.add(ref_chr);
                                                potentialinsertion.add(alignment_next.getRefEndC());
                                                potentialinsertion.add(alignment_current.getRefStartC());
                                                potentialinsertion.add(Math.abs(deviation));
                                                potentialinsertion.add("suppl");
                                                potentialinsertion.add("rev");
                                                potentialinsertion.add(confidence_level);
                                                potentialinsertions.add(potentialinsertion);

                                            }

                                        }

                                    }
                                    else if(deviation<=-parameter_setting.min_sv_size){

                                        if(distance_on_read<=parameter_setting.segment_gap_tolerance){

                                            if(alignment_current.getRefStrand().equals("+")){

                                                ArrayList potentialdeletion=new ArrayList();
                                                potentialdeletion.add(ref_chr);
                                                potentialdeletion.add(alignment_current.getRefEndC());
                                                potentialdeletion.add(alignment_next.getRefStartC());
                                                potentialdeletion.add(Math.abs(deviation));
                                                potentialdeletion.add("suppl");
                                                potentialdeletion.add("fwd");
                                                potentialdeletion.add(confidence_level);
                                                potentialdeletions.add(potentialdeletion);

                                                ArrayList potentialtranslocationbreakpoint=new ArrayList();
                                                potentialtranslocationbreakpoint.add("fwd");
                                                potentialtranslocationbreakpoint.add("fwd");
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefEndC());
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefStartC());
                                                potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                                potentialtranslocationbreakpoint.add(is_used);
                                                potentialtranslocationbreakpoint.add("del");
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add("del");
                                                potentialtranslocationbreakpoint.add(confidence_level);
                                                potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                            }
                                            else{

                                                ArrayList potentialdeletion=new ArrayList();
                                                potentialdeletion.add(ref_chr);
                                                potentialdeletion.add(alignment_next.getRefEndC());
                                                potentialdeletion.add(alignment_current.getRefStartC());
                                                potentialdeletion.add(Math.abs(deviation));
                                                potentialdeletion.add("suppl");
                                                potentialdeletion.add("rev");
                                                potentialdeletion.add(confidence_level);
                                                potentialdeletions.add(potentialdeletion);

                                                ArrayList potentialtranslocationbreakpoint=new ArrayList();
                                                potentialtranslocationbreakpoint.add("rev");
                                                potentialtranslocationbreakpoint.add("rev");
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefStartC());
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefEndC());
                                                potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                                potentialtranslocationbreakpoint.add(is_used);
                                                potentialtranslocationbreakpoint.add("del");
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add("del");
                                                potentialtranslocationbreakpoint.add(confidence_level);
                                                potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                            }

                                        }

                                    }

                                }
                                else{

                                    if(distance_on_reference<-parameter_setting.min_sv_size){

                                        if(alignment_current.getRefStrand().equals("+")){

                                            if(alignment_next.getRefEndC()>alignment_current.getRefStartC()){

                                                if(Math.abs(alignment_current.getScore()-alignment_next.getScore())/Math.max(alignment_current.getScore(),alignment_next.getScore())<0.5){

                                                    ArrayList tandem_duplication=new ArrayList();
                                                    tandem_duplication.add(ref_chr);
                                                    tandem_duplication.add(alignment_next.getRefStartC());
                                                    tandem_duplication.add(alignment_current.getRefEndC());
                                                    tandem_duplication.add(confidence_level);
                                                    tandem_duplications.add(tandem_duplication);

                                                }

                                            }
                                            else if(distance_on_reference>=-parameter_setting.max_sv_size){

                                                if(Math.abs(alignment_current.getScore()-alignment_next.getScore())/Math.max(alignment_current.getScore(),alignment_next.getScore())<0.5){

                                                    ArrayList tandem_duplication=new ArrayList();
                                                    tandem_duplication.add(ref_chr);
                                                    tandem_duplication.add(alignment_next.getRefStartC());
                                                    tandem_duplication.add(alignment_current.getRefEndC());
                                                    tandem_duplication.add(confidence_level);
                                                    tandem_duplications.add(tandem_duplication);

                                                }

                                                ArrayList potentialtranslocationbreakpoint =new ArrayList();
                                                potentialtranslocationbreakpoint.add("fwd");
                                                potentialtranslocationbreakpoint.add("fwd");
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefEndC());
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefStartC());
                                                potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                                potentialtranslocationbreakpoint.add(is_used);
                                                potentialtranslocationbreakpoint.add("intersect");
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add("intersect");
                                                potentialtranslocationbreakpoint.add(confidence_level);
                                                potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                            }
                                            else{

                                                ArrayList potentialtranslocationbreakpoint =new ArrayList();
                                                potentialtranslocationbreakpoint.add("fwd");
                                                potentialtranslocationbreakpoint.add("fwd");
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefEndC());
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefStartC());
                                                potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                                potentialtranslocationbreakpoint.add(is_used);
                                                potentialtranslocationbreakpoint.add("intersect");
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add("intersect");
                                                potentialtranslocationbreakpoint.add(confidence_level);
                                                potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                            }

                                        }
                                        else{

                                            if(alignment_next.getRefStartC()<alignment_current.getRefEndC()) {

                                                if(Math.abs(alignment_current.getScore()-alignment_next.getScore())/Math.max(alignment_current.getScore(),alignment_next.getScore())<0.5){

                                                    ArrayList tandem_duplication=new ArrayList();
                                                    tandem_duplication.add(ref_chr);
                                                    tandem_duplication.add(alignment_current.getRefStartC());
                                                    tandem_duplication.add(alignment_next.getRefEndC());
                                                    tandem_duplication.add(confidence_level);
                                                    tandem_duplications.add(tandem_duplication);

                                                }

                                            }
                                            else if(distance_on_reference>=-parameter_setting.max_sv_size){

                                                if(Math.abs(alignment_current.getScore()-alignment_next.getScore())/Math.max(alignment_current.getScore(),alignment_next.getScore())<0.5){

                                                    ArrayList tandem_duplication=new ArrayList();
                                                    tandem_duplication.add(ref_chr);
                                                    tandem_duplication.add(alignment_current.getRefStartC());
                                                    tandem_duplication.add(alignment_next.getRefEndC());
                                                    tandem_duplication.add(confidence_level);
                                                    tandem_duplications.add(tandem_duplication);

                                                }

                                                ArrayList potentialtranslocationbreakpoint=new ArrayList();
                                                potentialtranslocationbreakpoint.add("rev");
                                                potentialtranslocationbreakpoint.add("rev");
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefStartC());
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefEndC());
                                                potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                                potentialtranslocationbreakpoint.add(is_used);
                                                potentialtranslocationbreakpoint.add("intersect");
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add("intersect");
                                                potentialtranslocationbreakpoint.add(confidence_level);
                                                potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                            }
                                            else{

                                                ArrayList potentialtranslocationbreakpoint=new ArrayList();
                                                potentialtranslocationbreakpoint.add("rev");
                                                potentialtranslocationbreakpoint.add("rev");
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefStartC());
                                                potentialtranslocationbreakpoint.add(ref_chr);
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefEndC());
                                                potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                                potentialtranslocationbreakpoint.add(is_used);
                                                potentialtranslocationbreakpoint.add("intersect");
                                                potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                                potentialtranslocationbreakpoint.add("intersect");
                                                potentialtranslocationbreakpoint.add(confidence_level);
                                                potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                            }

                                        }

                                    }

                                }

                            }

                        }
                        else{

                            if(alignment_current.getRefStrand().equals("+")&&alignment_next.getRefStrand().equals("-")){

                                if(-parameter_setting.segment_overlap_tolerance<=distance_on_read&&distance_on_read<=parameter_setting.segment_gap_tolerance){

                                    if((alignment_next.getRefStartC()-alignment_current.getRefEndC())>=-parameter_setting.segment_overlap_tolerance){

                                        if((alignment_next.getRefEndC()-alignment_current.getRefEndC())<=parameter_setting.max_sv_size){

                                            Global.inversion_signatures.add(new SignatureInversion(ref_chr,alignment_current.getRefEndC(),alignment_next.getRefEndC(),"suppl",read_name,"left_fwd",confidence_level));

                                        }
                                        else{

                                            ArrayList potentialtranslocationbreakpoint=new ArrayList();
                                            potentialtranslocationbreakpoint.add("fwd");
                                            potentialtranslocationbreakpoint.add("rev");
                                            potentialtranslocationbreakpoint.add(ref_chr);
                                            potentialtranslocationbreakpoint.add(alignment_current.getRefEndC());
                                            potentialtranslocationbreakpoint.add(ref_chr);
                                            potentialtranslocationbreakpoint.add(alignment_next.getRefEndC());
                                            potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                            potentialtranslocationbreakpoint.add(is_used);
                                            potentialtranslocationbreakpoint.add("del");
                                            potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                            potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                            potentialtranslocationbreakpoint.add("inv");
                                            potentialtranslocationbreakpoint.add(confidence_level);
                                            potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                        }

                                    }
                                    else if(alignment_current.getRefStartC()-alignment_next.getRefEndC()>=-parameter_setting.segment_overlap_tolerance){

                                        if(alignment_current.getRefEndC()-alignment_next.getRefEndC()<=parameter_setting.max_sv_size){

                                            Global.inversion_signatures.add(new SignatureInversion(ref_chr,alignment_next.getRefEndC(),alignment_current.getRefEndC(),"suppl",read_name,"left_rev",confidence_level));

                                        }
                                        else{

                                            ArrayList potentialtranslocationbreakpoint=new ArrayList();
                                            potentialtranslocationbreakpoint.add("fwd");
                                            potentialtranslocationbreakpoint.add("rev");
                                            potentialtranslocationbreakpoint.add(ref_chr);
                                            potentialtranslocationbreakpoint.add(alignment_current.getRefEndC());
                                            potentialtranslocationbreakpoint.add(ref_chr);
                                            potentialtranslocationbreakpoint.add(alignment_next.getRefEndC());
                                            potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                            potentialtranslocationbreakpoint.add(is_used);
                                            potentialtranslocationbreakpoint.add("intersect");
                                            potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                            potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                            potentialtranslocationbreakpoint.add("inv");
                                            potentialtranslocationbreakpoint.add(confidence_level);
                                            potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                        }

                                    }

                                }

                            }

                            if(alignment_current.getRefStrand().equals("-")&&alignment_next.getRefStrand().equals("+")){

                                if(-parameter_setting.segment_overlap_tolerance<=distance_on_read&&distance_on_read<=parameter_setting.segment_gap_tolerance){

                                    if((alignment_next.getRefStartC()-alignment_current.getRefEndC())>=-parameter_setting.segment_overlap_tolerance){

                                        if((alignment_next.getRefStartC()-alignment_current.getRefStartC())<=parameter_setting.max_sv_size){

                                            Global.inversion_signatures.add(new SignatureInversion(ref_chr,alignment_current.getRefStartC(),alignment_next.getRefStartC(),"suppl",read_name,"right_fwd",confidence_level));

                                        }
                                        else{

                                            ArrayList potentialtranslocationbreakpoint=new ArrayList();
                                            potentialtranslocationbreakpoint.add("rev");
                                            potentialtranslocationbreakpoint.add("fwd");
                                            potentialtranslocationbreakpoint.add(ref_chr);
                                            potentialtranslocationbreakpoint.add(alignment_current.getRefStartC());
                                            potentialtranslocationbreakpoint.add(ref_chr);
                                            potentialtranslocationbreakpoint.add(alignment_next.getRefStartC());
                                            potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                            potentialtranslocationbreakpoint.add(is_used);
                                            potentialtranslocationbreakpoint.add("del");
                                            potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                            potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                            potentialtranslocationbreakpoint.add("inv");
                                            potentialtranslocationbreakpoint.add(confidence_level);
                                            potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                        }

                                    }
                                    else if(alignment_current.getRefStartC()-alignment_next.getRefEndC()>=-parameter_setting.segment_overlap_tolerance){

                                        if(alignment_current.getRefStartC()-alignment_next.getRefStartC()<=parameter_setting.max_sv_size){

                                            Global.inversion_signatures.add(new SignatureInversion(ref_chr,alignment_next.getRefStartC(),alignment_current.getRefStartC(),"suppl",read_name,"right_rev",confidence_level));

                                        }
                                        else{

                                            ArrayList potentialtranslocationbreakpoint=new ArrayList();
                                            potentialtranslocationbreakpoint.add("rev");
                                            potentialtranslocationbreakpoint.add("fwd");
                                            potentialtranslocationbreakpoint.add(ref_chr);
                                            potentialtranslocationbreakpoint.add(alignment_current.getRefStartC());
                                            potentialtranslocationbreakpoint.add(ref_chr);
                                            potentialtranslocationbreakpoint.add(alignment_next.getRefStartC());
                                            potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                            potentialtranslocationbreakpoint.add(is_used);
                                            potentialtranslocationbreakpoint.add("intersect");
                                            potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                            potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                            potentialtranslocationbreakpoint.add("inv");
                                            potentialtranslocationbreakpoint.add(confidence_level);
                                            potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                        }

                                    }

                                }

                            }

                        }

                    }
                    else{

                        String ref_chr_current= alignment_current.getRef();
                        String ref_chr_next=alignment_next.getRef();

                        if(alignment_current.getRefStrand().equals(alignment_next.getRefStrand())){

                            if(distance_on_read>=-parameter_setting.segment_overlap_tolerance){

                                if(distance_on_read<=parameter_setting.segment_gap_tolerance){

                                    if(alignment_current.getRefStrand().equals("+")){

                                        ArrayList potentialtranslocationbreakpoint=new ArrayList();
                                        potentialtranslocationbreakpoint.add("fwd");
                                        potentialtranslocationbreakpoint.add("fwd");
                                        potentialtranslocationbreakpoint.add(ref_chr_current);
                                        potentialtranslocationbreakpoint.add(alignment_current.getRefEndC());
                                        potentialtranslocationbreakpoint.add(ref_chr_next);
                                        potentialtranslocationbreakpoint.add(alignment_next.getRefStartC());
                                        potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                        potentialtranslocationbreakpoint.add(is_used);
                                        potentialtranslocationbreakpoint.add("intertrans");
                                        potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                        potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                        potentialtranslocationbreakpoint.add("intertrans");
                                        potentialtranslocationbreakpoint.add(confidence_level);
                                        potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                    }
                                    else{
                                        ArrayList potentialtranslocationbreakpoint=new ArrayList();
                                        potentialtranslocationbreakpoint.add("rev");
                                        potentialtranslocationbreakpoint.add("rev");
                                        potentialtranslocationbreakpoint.add(ref_chr_current);
                                        potentialtranslocationbreakpoint.add(alignment_current.getRefStartC());
                                        potentialtranslocationbreakpoint.add(ref_chr_next);
                                        potentialtranslocationbreakpoint.add(alignment_next.getRefEndC());
                                        potentialtranslocationbreakpoint.add(alignment_current.getReadStartC());
                                        potentialtranslocationbreakpoint.add(is_used);
                                        potentialtranslocationbreakpoint.add("intertrans");
                                        potentialtranslocationbreakpoint.add(alignment_current.getRefAlignLen());
                                        potentialtranslocationbreakpoint.add(alignment_next.getRefAlignLen());
                                        potentialtranslocationbreakpoint.add("intertrans");
                                        potentialtranslocationbreakpoint.add(confidence_level);
                                        potentialtranslocationsbreakpoint.add(potentialtranslocationbreakpoint);

                                    }

                                }

                            }

                        }

                    }

                }

                String current_chromosome=null;
                List current_starts=new LinkedList();
                List current_ends=new LinkedList();
                List current_confidence_level=new LinkedList();
                int currrent_copy_number=0;

                if(tandem_duplications!=null){

                    int tandem_duplications_size=tandem_duplications.size();

                    for(int n=0;n<tandem_duplications_size;n++){

                        ArrayList tandem_duplication=tandem_duplications.get(n);

                        if(current_chromosome==null){

                            current_chromosome=(String)(tandem_duplication.get(0));
                            current_starts.add(tandem_duplication.get(1));
                            current_ends.add(tandem_duplication.get(2));
                            current_confidence_level.add(tandem_duplication.get(3));
                            currrent_copy_number=1;
                        }
                        else{

                            if(is_similar(current_chromosome,(long)mean(current_starts),(long)mean(current_ends),(String)tandem_duplication.get(0),(long)tandem_duplication.get(1),(long)tandem_duplication.get(2))){

                                current_starts.add(tandem_duplication.get(1));
                                current_ends.add(tandem_duplication.get(2));
                                current_confidence_level.add(tandem_duplication.get(3));
                                currrent_copy_number+=1;

                            }
                            else{

                                Global.duplication_tandem_signatures.add(new SignatureDuplicationTandem(current_chromosome,(long)mean(current_starts),(long)mean(current_ends),currrent_copy_number,"suppl",read_name,(double)average(current_confidence_level)));
                                current_chromosome=(String) tandem_duplication.get(0);
                                current_starts.add(tandem_duplication.get(1));
                                current_ends.add(tandem_duplication.get(2));
                                current_confidence_level.add(tandem_duplication.get(3));
                                currrent_copy_number=1;

                            }

                        }

                    }

                }

                if(current_chromosome!=null){

                    Global.duplication_tandem_signatures.add(new SignatureDuplicationTandem(current_chromosome,(long)mean(current_starts),(long)mean(current_ends),currrent_copy_number,"suppl",read_name,(double)average(current_confidence_level)));

                }

                if(potentialtranslocationsbreakpoint!=null){

                    Collections.sort(potentialtranslocationsbreakpoint,new potentialtranslocationsbreakpointreadCorOrder());
                    int potentialtranslocationsbreakpoint_size=potentialtranslocationsbreakpoint.size();

                    if(potentialtranslocationsbreakpoint_size>0){

                        for(int n=0;n<potentialtranslocationsbreakpoint_size-1;n++){

                            ArrayList potentialtranslocationsbreakpoint_current=potentialtranslocationsbreakpoint.get(n);
                            String dir1_current=(String) potentialtranslocationsbreakpoint_current.get(0);
                            String dir2_current=(String) potentialtranslocationsbreakpoint_current.get(1);
                            String chr1_current=(String) potentialtranslocationsbreakpoint_current.get(2);
                            long pos1_current=(long) potentialtranslocationsbreakpoint_current.get(3);
                            String chr2_current=(String) potentialtranslocationsbreakpoint_current.get(4);
                            long pos2_current=(long) potentialtranslocationsbreakpoint_current.get(5);
                            String type_current=(String) potentialtranslocationsbreakpoint_current.get(8);
                            long length1_current=(long)potentialtranslocationsbreakpoint_current.get(9);
                            long length2_current=(long)potentialtranslocationsbreakpoint_current.get(10);

                            ArrayList potentialtranslocationsbreakpoint_next=potentialtranslocationsbreakpoint.get(n+1);
                            String dir1_next=(String) potentialtranslocationsbreakpoint_next.get(0);
                            String dir2_next=(String) potentialtranslocationsbreakpoint_next.get(1);
                            String chr1_next=(String) potentialtranslocationsbreakpoint_next.get(2);
                            long pos1_next=(long) potentialtranslocationsbreakpoint_next.get(3);
                            String chr2_next=(String) potentialtranslocationsbreakpoint_next.get(4);
                            long pos2_next=(long) potentialtranslocationsbreakpoint_next.get(5);
                            String type_next=(String) potentialtranslocationsbreakpoint_next.get(8);
                            long length1_next=(long)potentialtranslocationsbreakpoint_next.get(9);
                            long length2_next=(long)potentialtranslocationsbreakpoint_next.get(10);

                            double confidence_level=((double)potentialtranslocationsbreakpoint_current.get(12)+(double)potentialtranslocationsbreakpoint_next.get(12))/2;

                            if(parameter_setting.min_sv_size<=Math.abs((pos2_current-pos1_next))&&Math.abs((pos2_current-pos1_next))<=parameter_setting.max_sv_size){

                                if(chr2_current.equals(chr1_next)&&dir2_current.equals(dir1_next)){

                                    if(chr1_current.equals(chr2_next)){

                                        if(chr1_current.equals(chr1_next)){

                                            if(!type_current.equals(type_next)){

                                                if(dir1_current.equals(dir2_current)&&dir2_next.equals(dir2_current)){

                                                    if(dir2_current.equals("fwd")){

                                                        Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos2_current,pos1_next,chr1_current,pos1_current,pos2_next,"suppl",read_name,confidence_level));
                                                        potentialtranslocationsbreakpoint_current.set(7,true);
                                                        potentialtranslocationsbreakpoint_next.set(7,true);

                                                    }
                                                    else if(dir2_current.equals("rev")){

                                                        Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos1_next,pos2_current,chr1_current,pos2_next,pos1_current,"suppl",read_name,confidence_level));
                                                        potentialtranslocationsbreakpoint_current.set(7,true);
                                                        potentialtranslocationsbreakpoint_next.set(7,true);

                                                    }

                                                }
                                                else if(dir1_current.equals(dir2_current)&&!dir2_next.equals(dir2_current)){

                                                    if(dir2_current.equals("fwd")){

                                                        Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos2_current,pos1_next,chr1_current,pos1_current,pos2_next-length2_next,"suppl",read_name,confidence_level));
                                                        potentialtranslocationsbreakpoint_current.set(7,true);
                                                        potentialtranslocationsbreakpoint_next.set(7,true);

                                                    }
                                                    else if(dir2_current.equals("rev")){

                                                        Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos1_next,pos2_current,chr1_current,pos2_next+length2_next,pos1_current,"suppl",read_name,confidence_level));
                                                        potentialtranslocationsbreakpoint_current.set(7,true);
                                                        potentialtranslocationsbreakpoint_next.set(7,true);

                                                    }

                                                }
                                                else if(!dir1_current.equals(dir2_current)&&dir2_next.equals(dir2_current)){

                                                    if(dir2_current.equals("fwd")){

                                                        Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos2_current,pos1_next,chr1_current,pos1_current+length1_current,pos2_next,"suppl",read_name,confidence_level));
                                                        potentialtranslocationsbreakpoint_current.set(7,true);
                                                        potentialtranslocationsbreakpoint_next.set(7,true);

                                                    }
                                                    else if(dir2_current.equals("rev")){

                                                        Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos1_next,pos2_current,chr1_current,pos2_next,pos1_current-length1_current,"suppl",read_name,confidence_level));
                                                        potentialtranslocationsbreakpoint_current.set(7,true);
                                                        potentialtranslocationsbreakpoint_next.set(7,true);

                                                    }

                                                }
                                                else if(!dir1_current.equals(dir2_current)&&!dir2_next.equals(dir2_current)){

                                                    if(dir2_current.equals("fwd")){

                                                        Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos2_current,pos1_next,chr1_current,pos2_next,pos1_current,"suppl",read_name,confidence_level));
                                                        potentialtranslocationsbreakpoint_current.set(7,true);
                                                        potentialtranslocationsbreakpoint_next.set(7,true);

                                                    }
                                                    else if(dir2_current.equals("rev")){

                                                        Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos1_next,pos2_current,chr1_current,pos1_current,pos2_next,"suppl",read_name,confidence_level));
                                                        potentialtranslocationsbreakpoint_current.set(7,true);
                                                        potentialtranslocationsbreakpoint_next.set(7,true);

                                                    }

                                                }

                                            }

                                        }
                                        else{

                                            if(dir1_current.equals(dir2_current)&&dir2_next.equals(dir2_current)){

                                                if(dir2_current.equals("fwd")){

                                                    Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos2_current,pos1_next,chr1_current,pos1_current,pos2_next,"suppl",read_name,confidence_level));
                                                    potentialtranslocationsbreakpoint_current.set(7,true);
                                                    potentialtranslocationsbreakpoint_next.set(7,true);

                                                }
                                                else if(dir2_current.equals("rev")){

                                                    Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos1_next,pos2_current,chr1_current,pos2_next,pos1_current,"suppl",read_name,confidence_level));
                                                    potentialtranslocationsbreakpoint_current.set(7,true);
                                                    potentialtranslocationsbreakpoint_next.set(7,true);

                                                }

                                            }
                                            else if(dir1_current.equals(dir2_current)&&!dir2_next.equals(dir2_current)){

                                                if(dir2_current.equals("fwd")){

                                                    Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos2_current,pos1_next,chr1_current,pos1_current,pos2_next-length2_next,"suppl",read_name,confidence_level));
                                                    potentialtranslocationsbreakpoint_current.set(7,true);
                                                    potentialtranslocationsbreakpoint_next.set(7,true);

                                                }
                                                else if(dir2_current.equals("rev")){

                                                    Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos1_next,pos2_current,chr1_current,pos2_next+length2_next,pos1_current,"suppl",read_name,confidence_level));
                                                    potentialtranslocationsbreakpoint_current.set(7,true);
                                                    potentialtranslocationsbreakpoint_next.set(7,true);

                                                }

                                            }
                                            else if(!dir1_current.equals(dir2_current)&&dir2_next.equals(dir2_current)){

                                                if(dir2_current.equals("fwd")){

                                                    Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos2_current,pos1_next,chr1_current,pos1_current+length1_current,pos2_next,"suppl",read_name,confidence_level));
                                                    potentialtranslocationsbreakpoint_current.set(7,true);
                                                    potentialtranslocationsbreakpoint_next.set(7,true);

                                                }
                                                else if(dir2_current.equals("rev")){

                                                    Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos1_next,pos2_current,chr1_current,pos2_next,pos1_current-length1_current,"suppl",read_name,confidence_level));
                                                    potentialtranslocationsbreakpoint_current.set(7,true);
                                                    potentialtranslocationsbreakpoint_next.set(7,true);

                                                }

                                            }
                                            else if(!dir1_current.equals(dir2_current)&&!dir2_next.equals(dir2_current)){

                                                if(dir2_current.equals("fwd")){

                                                    Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos2_current,pos1_next,chr1_current,pos2_next,pos1_current,"suppl",read_name,confidence_level));
                                                    potentialtranslocationsbreakpoint_current.set(7,true);
                                                    potentialtranslocationsbreakpoint_next.set(7,true);

                                                }
                                                else if(dir2_current.equals("rev")){

                                                    Global.Potential_Trans_Inters_signatures.add(new SignaturePotential_Trans_Inters(chr2_current,pos1_next,pos2_current,chr1_current,pos1_current,pos2_next,"suppl",read_name,confidence_level));
                                                    potentialtranslocationsbreakpoint_current.set(7,true);
                                                    potentialtranslocationsbreakpoint_next.set(7,true);

                                                }

                                            }

                                        }

                                    }

                                }

                            }

                        }

                        for(int i=0;i<potentialtranslocationsbreakpoint_size;i++){

                            ArrayList potentialtranslocationbreakpoint =potentialtranslocationsbreakpoint.get(i);

                            if(potentialtranslocationbreakpoint.get(11).equals("del")){

                                if((boolean)potentialtranslocationbreakpoint.get(7)){

                                    int potentialdeletions_size=potentialdeletions.size();

                                    if(potentialdeletions_size>0){

                                        for(int j=0;j<potentialdeletions.size();j++){

                                            ArrayList potentialdeletion=potentialdeletions.get(j);

                                            if(potentialdeletion.get(5).equals("fwd")){

                                                if(potentialtranslocationbreakpoint.get(2).equals(potentialdeletion.get(0))&&(long)potentialtranslocationbreakpoint.get(3)==(long)potentialdeletion.get(1)&&(long)potentialtranslocationbreakpoint.get(5)==(long)potentialdeletion.get(2)) {

                                                    potentialdeletions.remove(j);
                                                    j--;

                                                }

                                            }
                                            else if(potentialdeletion.get(5).equals("rev")){

                                                if(potentialtranslocationbreakpoint.get(2).equals(potentialdeletion.get(0))&&(long)potentialtranslocationbreakpoint.get(3)==(long)potentialdeletion.get(2)&&(long)potentialtranslocationbreakpoint.get(5)==(long)potentialdeletion.get(1)) {

                                                    potentialdeletions.remove(j);
                                                    j--;

                                                }

                                            }

                                        }

                                    }

                                }

                            }

                            if(potentialtranslocationbreakpoint.get(11).equals("intersect")){

                                if((boolean)potentialtranslocationbreakpoint.get(7)){

                                    int tandemduplication_size=Global.duplication_tandem_signatures.size();

                                    if(tandemduplication_size>0){

                                        for(int j=0;j<Global.duplication_tandem_signatures.size();j++){

                                            SignatureDuplicationTandem signatureDuplicationTandem=Global.duplication_tandem_signatures.get(j);

                                            if(potentialtranslocationbreakpoint.get(2).equals(signatureDuplicationTandem.contig)&&(long)potentialtranslocationbreakpoint.get(3)==signatureDuplicationTandem.start_position_on_ref&&(long)potentialtranslocationbreakpoint.get(5)==signatureDuplicationTandem.end_position_on_ref) {

                                                Global.duplication_tandem_signatures.remove(j);
                                                j--;

                                            }
                                            else if(potentialtranslocationbreakpoint.get(2).equals(signatureDuplicationTandem.contig)&&(long)potentialtranslocationbreakpoint.get(3)==signatureDuplicationTandem.end_position_on_ref&&(long)potentialtranslocationbreakpoint.get(5)==signatureDuplicationTandem.start_position_on_ref){

                                                Global.duplication_tandem_signatures.remove(j);
                                                j--;

                                            }

                                        }

                                    }

                                }

                            }

                            Global.TranslocationBreakpoint_signatures.add(new SignatureTranslocationBreakpoint((String) potentialtranslocationbreakpoint.get(2),(long)potentialtranslocationbreakpoint.get(3),(String) potentialtranslocationbreakpoint.get(0),(String) potentialtranslocationbreakpoint.get(4),(long)potentialtranslocationbreakpoint.get(5),(String) potentialtranslocationbreakpoint.get(1),"suppl",(double)potentialtranslocationbreakpoint.get(12)));

                        }

                    }

                }

                if(potentialinsertions!=null){

                    int size=potentialinsertions.size();
                    for(int i=0;i<size;i++){

                        ArrayList potentialinsertion=potentialinsertions.get(i);
                        Global.insertion_signatures.add(new SignatureInsertion((String)potentialinsertion.get(0),(long)potentialinsertion.get(1),(long)potentialinsertion.get(2),(String) potentialinsertion.get(4),(String) potentialinsertion.get(5),read_name,(double)potentialinsertion.get(6),(long)potentialinsertion.get(3)));

                    }

                }

                if(potentialdeletions!=null){

                    int size=potentialdeletions.size();
                    for(int i=0;i<size;i++){

                        ArrayList potentialdeletion=potentialdeletions.get(i);
                        Global.deletion_signatures.add(new SignatureDeletion((String)potentialdeletion.get(0),(long)potentialdeletion.get(1),(long)potentialdeletion.get(2),(String) potentialdeletion.get(4),(String) potentialdeletion.get(5),read_name,(double)potentialdeletion.get(6),(long)potentialdeletion.get(3)));

                    }

                }

            }

        }

    }

    private static Object average(List current_confidence_level) {

        double current_position_sum=0;

        for(int n=0;n<current_confidence_level.size();n++){

            current_position_sum+=(double)current_confidence_level.get(n);

        }

        return current_position_sum/current_confidence_level.size();

    }

    private static boolean is_similar(String chr1, long start1, long end1, String chr2, long start2, long end2) {

        if(chr1.equals(chr2)&&Math.abs(start1-start2)<20&&Math.abs(end1-end2)<20){

            return true;

        }
        else{

            return false;

        }

    }

    private static Object mean(List current_positions) {

        long current_position_sum=0;

        for(int n=0;n<current_positions.size();n++){

            current_position_sum+=(long)current_positions.get(n);

        }

        return current_position_sum/current_positions.size();

    }

}