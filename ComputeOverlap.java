

import java.util.Arrays;

public class ComputeOverlap {

    private long overlap;
    private long seedlen;
    private long alignlen;

    private void setOverlap(long overlap){
        this.overlap=overlap;
    }

    public long getOverlap(){
        return overlap;
    }

    private void setSeedlen(long seedlen){
        this.seedlen=seedlen;
    }

    public long getSeedlen(){
        return seedlen;
    }

    private void setAlignlen(long alignlen){
        this.alignlen=alignlen;
    }

    public long getAlignlen(){
        return alignlen;
    }

    public static ComputeOverlap computeoverlap(long Segment1S,long Segment1E,long Segment2S,long Segment2E){

        long order[]={Segment1S,Segment1E,Segment2S,Segment2E};
        ComputeOverlap result=new ComputeOverlap();

        if(Segment2S>Segment1E){
            result.setOverlap(0);
            result.setSeedlen(Segment1E-Segment1S+1);
            result.setAlignlen(Segment2E-Segment2S+1);
        }
        else if(Segment1S>Segment2E){
            result.setOverlap(0);
            result.setSeedlen(Segment1E-Segment1S+1);
            result.setAlignlen(Segment2E-Segment2S+1);
        }
        else {
            Arrays.sort(order);
            result.setOverlap((Segment1E-Segment1S+1)+(Segment2E-Segment2S+1)-(order[3]-order[0]+1));
            result.setSeedlen(Segment1E-Segment1S+1);
            result.setAlignlen(Segment2E-Segment2S+1);
        }
        return result;
    }
}
