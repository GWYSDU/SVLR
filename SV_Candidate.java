
import java.util.*;
public class SV_Candidate {

    String contig;
    long start_position_on_ref;
    long end_position_on_ref;
    double score;

    public SV_Candidate(){

    }

}

class SV_CandidateDeletion extends SV_Candidate{

    double span;

    public SV_CandidateDeletion(String contig, long start_position_on_ref, long end_position_on_ref,double score,double span) {
        super();
        this.contig =contig;
        this.start_position_on_ref = start_position_on_ref;
        this.end_position_on_ref = end_position_on_ref;
        this.score=score;
        this.span=span;

    }

}

class SV_CandidateInsertion extends SV_Candidate{

    double span;
    public SV_CandidateInsertion(String contig, long start_position_on_ref, long end_position_on_ref,double score,double span) {
        super();
        this.contig =contig;
        this.start_position_on_ref = start_position_on_ref;
        this.end_position_on_ref = end_position_on_ref;
        this.score=score;
        this.span=span;

    }
}

class SV_CandidateInversion extends SV_Candidate{

    public SV_CandidateInversion(String contig, long start_position_on_ref, long end_position_on_ref,double score) {
        super();
        this.contig =contig;
        this.start_position_on_ref = start_position_on_ref;
        this.end_position_on_ref = end_position_on_ref;
        this.score=score;

    }

}

class SV_CandidateDuplicationTandem extends SV_Candidate{

    int copies;

    public SV_CandidateDuplicationTandem(String contig, long start, long end,int copies,double score) {

        super();
        this.contig =contig;
        this.start_position_on_ref = start;
        this.end_position_on_ref = end;
        this.copies=copies;
        this.score=score;

    }

}

class SV_CandidateTranslocation extends SV_Candidate{

    String contig1;
    String contig2;
    long destination_start_position_on_ref;
    long destination_end_position_on_ref;

    public SV_CandidateTranslocation( String contig1,long start_position_on_ref,long end_position_on_ref,String contig2,long destination_start_position_on_ref,long destination_end_position_on_ref,double score){
        super();
        this.start_position_on_ref = start_position_on_ref;
        this.end_position_on_ref = end_position_on_ref;
        this.contig1=contig1;
        this.contig2=contig2;
        this.destination_start_position_on_ref=destination_start_position_on_ref;
        this.destination_end_position_on_ref=destination_end_position_on_ref;
        this.score=score;

    }

}
class SV_CandidateInterspersedDuplication extends SV_Candidate{

    String contig1;
    String contig2;
    long destination_start_position_on_ref;
    long destination_end_position_on_ref;
    boolean is_cutpaste;

    public SV_CandidateInterspersedDuplication( String contig1,long start_position_on_ref,long end_position_on_ref,String contig2,long destination_start_position_on_ref,long destination_end_position_on_ref,double score,boolean is_cutpaste){

        super();
        this.start_position_on_ref = start_position_on_ref;
        this.end_position_on_ref = end_position_on_ref;
        this.contig1=contig1;
        this.contig2=contig2;
        this.destination_start_position_on_ref=destination_start_position_on_ref;
        this.destination_end_position_on_ref=destination_end_position_on_ref;
        this.score=score;
        this.is_cutpaste=is_cutpaste;

    }

}
