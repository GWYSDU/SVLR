
import java.util.*;
public class SVsignature {

    long span;
    String contig;
    long start_position_on_ref;
    long end_position_on_ref;
    long start_position_on_read;
    String signature;
    String direction;
    String read;
    String type;
    double confidence_level;

    public SVsignature(){

    }

    public String getSourceContig(){

        return contig;

    }

    public long getSourceStart(){

        return start_position_on_ref;

    }
    public long getSourceEnd(){

        return end_position_on_ref;

    }
    public String getType(){

        return type;

    }

    public String get_key_type(){

        return getType();

    }
    public String get_key_contig(){

        return getSourceContig();

    }
    public long  get_key_position(){

        return (getSourceStart()+getSourceEnd())/2;

    }
    public long getSourceMidpoint(){

        return (getSourceStart()+getSourceEnd())/2;

    }

    public double getConfidence_level(){

        return confidence_level;

    }

    public double mean_distance_to(SVsignature SVsignature2){

        if(this.getType().equals(SVsignature2.getType())&&this.getSourceContig().equals(SVsignature2.getSourceContig())){

            return Math.abs((this.getSourceMidpoint())-(SVsignature2.getSourceMidpoint()));

        }
        else{

            return Double.POSITIVE_INFINITY;

        }

    }

    public double span_loc_distance(SVsignature SVsignature2,double distance_normalizer){


        if(!this.getSourceContig().equals(SVsignature2.getSourceContig())){

            return Double.POSITIVE_INFINITY;

        }

        long this_span=this.getSourceEnd()-this.getSourceStart();
        long other_span=SVsignature2.getSourceEnd()-SVsignature2.getSourceStart();

        if(this.type.equals("ins")||this.type.equals("del")){

            this_span=this.span;

        }
        if(SVsignature2.type.equals("ins")||SVsignature2.type.equals("del")){

            other_span=SVsignature2.span;

        }

        double dist_span=Math.abs(this_span-other_span)/Math.max(this_span,other_span);
        long this_center=this.get_key_position();
        long other_center=SVsignature2.get_key_position();
        double dist_loc=Math.min(Math.min(Math.abs(this.getSourceStart()-SVsignature2.getSourceStart()),Math.abs(this.getSourceEnd()-SVsignature2.getSourceEnd())),Math.abs(this_center-other_center)/distance_normalizer);
        return dist_span+dist_loc;

    }

}

class SignatureDeletion extends SVsignature{

    public SignatureDeletion(String contig, long start_position_on_ref, long end_position_on_ref,String signature,String direction, String read, double confidence_level,long span) {
        super();
        this.contig =contig;
        this.start_position_on_ref = start_position_on_ref;
        this.end_position_on_ref = end_position_on_ref;
        this.signature=signature;
        this.direction=direction;
        this.read=read;
        this.confidence_level=confidence_level;
        this.type="del";
        this.span=span;

    }

}

class SignatureInsertion extends SVsignature{

    public SignatureInsertion(String contig, long start_position_on_ref, long end_position_on_ref,String signature,String direction, String read, double confidence_level,long span) {
        super();
        this.contig =contig;
        this.start_position_on_ref = start_position_on_ref;
        this.end_position_on_ref = end_position_on_ref;
        this.signature=signature;
        this.direction=direction;
        this.read=read;
        this.confidence_level=confidence_level;
        this.type="ins";
        this.span=span;
    }
}

class SignatureInversion extends SVsignature{

    public SignatureInversion(String contig, long start_position_on_ref, long end_position_on_ref,String signature,String read ,String direction,double cnfidence_level) {
        super();
        this.contig =contig;
        this.start_position_on_ref = start_position_on_ref;
        this.end_position_on_ref = end_position_on_ref;
        this.signature=signature;
        this.read=read;
        this.confidence_level=cnfidence_level;
        this.type="inv";
        this.direction=direction;
    }

}

class SignatureDuplicationTandem extends SVsignature{

    int copies;

    public SignatureDuplicationTandem(String contig, long start, long end,int copies, String signature, String read, double confidence_level) {

        super();
        this.contig =contig;
        this.start_position_on_ref = start;
        this.end_position_on_ref = end;
        this.signature=signature;
        this.read=read;
        this.confidence_level=confidence_level;
        this.copies=copies;
        this.type="dup";

    }

    public String getDestinationContig(){

        return getSourceContig();

    }

    public long getDestinationStart(){

        return getSourceEnd();

    }
    public long getDestinationEnd(){

        return getSourceEnd()+copies*(getSourceEnd()-getSourceStart());

    }

}

class SignatureTranslocationBreakpoint extends SVsignature{

    String contig1;
    long pos1;
    String direction1;
    String contig2;
    long pos2;
    String direction2;

    public SignatureTranslocationBreakpoint(String contig1,long pos1,String direction1,String contig2,long pos2,String direction2,String signature,double confidence_level) {

        super();
        this.signature=signature;
        this.confidence_level=confidence_level;
        this.contig1=contig1;
        this.pos1=pos1;
        this.direction1=direction1;
        this.contig2=contig2;
        this.pos2=pos2;
        this.direction2=direction2;
        this.type="trabp";

    }

    public String getSourceContig(){

        return contig1;

    }

    public long getSourceStart(){

        return pos1;

    }

    public long getSourceEnd(){

        return pos1+1;

    }

    public String getDestinationContig(){

        return contig2;

    }

    public long getDestinationStart(){

        return pos2;

    }

    public long getDestinationEnd(){

        return pos2+1;

    }

    public String get_key_contig(){

        return getSourceContig();

    }

    public long  get_key_position(){

        return getSourceStart();

    }

    public double mean_distance_to(SignatureTranslocationBreakpoint SVsignature2){


        if(this.getType().equals(SVsignature2.getType())&&this.getSourceContig().equals(SVsignature2.getSourceContig())){

            return Math.abs(this.get_key_position()-SVsignature2.get_key_position());

        }
        else{

            return Double.POSITIVE_INFINITY;

        }

    }

}

class SignaturePotential_Trans_Inters extends SVsignature{

    String contig1;
    String contig2;
    long destination_start_position_on_ref;
    long destination_end_position_on_ref;

    public SignaturePotential_Trans_Inters( String contig1,long start_position_on_ref,long end_position_on_ref,String contig2,long destination_start_position_on_ref,long destination_end_position_on_ref,String signature,String read,double confidence_level){

        super();
        this.start_position_on_ref = start_position_on_ref;
        this.end_position_on_ref = end_position_on_ref;
        this.signature=signature;
        this.read=read;
        this.confidence_level=confidence_level;
        this.contig1=contig1;
        this.contig2=contig2;
        this.destination_start_position_on_ref=destination_start_position_on_ref;
        this.destination_end_position_on_ref=destination_end_position_on_ref;
        this.type="potential_int_trans";

    }

    public String getSourceContig(){

        return contig1;

    }
    public long getSourceStart(){

        return start_position_on_ref;

    }
    public long getSourceEnd(){

        return end_position_on_ref;

    }
    public String getDestinationContig(){

        return contig2;

    }

    public long getDestinationStart(){

        return destination_start_position_on_ref;

    }
    public long getDestinationEnd(){

        return destination_end_position_on_ref;

    }

    public String get_key_contig(){

        return getSourceContig();

    }

    public String get_key_dest_contig(){

        return getDestinationContig();

    }

    public long  get_key_position(){

//        return getDestinationStart()+getSourceMidpoint();
        return getDestinationMidpoint();
    }

    public long getSourceMidpoint(){

        return (getSourceStart()+getSourceEnd())/2;

    }
    public long getDestinationMidpoint(){

        return (getDestinationStart()+getDestinationEnd())/2;

    }
    public double mean_distance_to(SignaturePotential_Trans_Inters SVsignature2){


        if(this.getSourceContig().equals(SVsignature2.getSourceContig())&&this.getDestinationContig().equals(SVsignature2.getDestinationContig())){

            return Math.abs(this.getSourceMidpoint()-SVsignature2.getSourceMidpoint())+Math.abs(this.getDestinationMidpoint()-SVsignature2.getDestinationMidpoint());

        }
        else{

            return Double.POSITIVE_INFINITY;

        }

    }

}

class SignatureClusterUniLocal extends  SVsignature{

    double score;
    int size;
    List<?>  members;
    double std_span;
    double std_pos;
    double ave_span;

    public SignatureClusterUniLocal(String contig, long start, long end, double score,int size, List<?> members,String type,double std_span,
                                    double std_pos,double ave_span) {

        super();
        this.contig =contig;
        this.start_position_on_ref = start;
        this.end_position_on_ref = end;
        this.score=score;
        this.std_span=std_span;
        this.std_pos=std_pos;
        this.size=size;
        this.members =members;
        this.type=type;
        this.ave_span=ave_span;

    }

    public long get_length(){

        return getSourceEnd()-getSourceStart();

    }

    public boolean equal(SignatureClusterUniLocal signatureClusterUniLocal) {

        if( this.contig.equals(signatureClusterUniLocal.contig)&&this.start_position_on_ref ==signatureClusterUniLocal.start_position_on_ref &&this.end_position_on_ref == signatureClusterUniLocal.end_position_on_ref &&this.score==signatureClusterUniLocal.score&&this.std_span==signatureClusterUniLocal.std_span&&this.std_pos==signatureClusterUniLocal.std_pos&&this.size==signatureClusterUniLocal.size&&this.members ==signatureClusterUniLocal.members&&this.type.equals(signatureClusterUniLocal.type)){

            return true;

        }
        else{

            return false;

        }

    }
}

class SignatureClusterBilocal extends SVsignature{

    String source_contig;
    long source_start;
    long source_end;
    String dest_contig;
    long dest_start;
    long dest_end;
    double score;
    int size;
    List<?> members;
    Double std_span;
    Double std_pos;

    public SignatureClusterBilocal(String source_contig,long source_start,long source_end,String dest_contig,long dest_start,long dest_end,double score,int size, List<?> members,String type,Double std_span,
                                   Double std_pos) {

        super();
        this.source_contig=source_contig;
        this.source_start=source_start;
        this.source_end=source_end;
        this.dest_contig=dest_contig;
        this.dest_start=dest_start;
        this.dest_end=dest_end;
        this.score=score;
        this.std_span=std_span;
        this.std_pos=std_pos;
        this.size=size;
        this.members=members;
        this.type=type;

    }

    public String getSourceContig(){

        return source_contig;

    }

    public long getSourceStart(){

        return source_start;

    }

    public long getSourceEnd(){

        return source_end;

    }

    public String getDestinationContig(){

        return dest_contig;

    }

    public long getDestinationStart(){

        return dest_start;

    }

    public long getDestinationEnd(){

        return dest_end;

    }

    public long get_source_length(){

        return getSourceEnd()-getSourceStart();

    }

    public long get_destination_length(){

        return getDestinationEnd()-getDestinationStart();

    }

}
