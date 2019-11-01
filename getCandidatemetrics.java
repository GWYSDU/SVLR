
import java.util.Arrays;
import java.util.List;

public class getCandidatemetrics {

    private long cscore;
    private double coveragefratction;
    private double per_base_score;

    private void setCscore(long cscore){
        this.cscore=cscore;
    }

    public long getCscore(){

        return cscore;
    }

    private void setCoveragefratction(double coveragefratction){
        this.coveragefratction=coveragefratction;
    }

    public double getCoveragefratction(){
        return coveragefratction;
    }
    public void setPer_base_score(double per_base_score){

        this.per_base_score=per_base_score;

    }

    public double getPer_base_score(){

        return per_base_score;

    }

    private long bestcscore;
    private double bestcoveragefratction;

    private void setbestCscore(long bestcscore){
        this.bestcscore=bestcscore;
    }

    public long getbestCscore(){

        return cscore;
    }

    private void setbestCoveragefratction(double bestcoveragefratction){
        this.bestcoveragefratction=bestcoveragefratction;
    }

    public double getbestCoveragefratction(){
        return bestcoveragefratction;
    }

    public static getCandidatemetrics getCandidatemetrics(List<parse_mapped_data> CandidateExtender){

        getCandidatemetrics result=new getCandidatemetrics();

        int CandidateExtenderNum = CandidateExtender.size();
        double coverageFraction = 0;
        double coverage=0;
        long startpos;
        long endpos;
        long cscore=0;

//        for (int i = 0; i < CandidateExtenderNum; i++){
//
//            cscore+=CandidateExtender.get(i).getScore();
//
//        }

        if (CandidateExtenderNum == 0) {

        } else if (CandidateExtenderNum == 1) {

            coverage= CandidateExtender.get(0).getReadAlignLen();
            cscore=CandidateExtender.get(0).getScore();
            startpos = CandidateExtender.get(0).getReadStartC();
            endpos = CandidateExtender.get(0).getReadEndC();

        } else {

            long x1s = CandidateExtender.get(0).getReadStartC();
            long x1e = CandidateExtender.get(0).getReadEndC();
            cscore=CandidateExtender.get(0).getScore();

            coverage= x1e - x1s + 1;

            for (int i = 1; i < CandidateExtenderNum; i++) {

                parse_mapped_data candidateRef = CandidateExtender.get(i);
                long overlap;
                long seedlen;
                long alignlen;
                long reduced_Score;

                ComputeOverlap result1;
                result1 = ComputeOverlap.computeoverlap(x1s, x1e, candidateRef.getReadStartC(), candidateRef.getReadEndC());
                overlap = result1.getOverlap();
                seedlen = result1.getSeedlen();
                alignlen = result1.getAlignlen();

                if (0 == overlap) {

                    coverage += alignlen;
                    cscore+=candidateRef.getScore();
//                    long gap=(x1e<candidateRef.getReadStartC())?(candidateRef.getReadStartC()-x1e):(x1s-candidateRef.getReadEndC());
//                    cscore-=gap;

                } else {

                    coverage += (alignlen - overlap);
                    reduced_Score=ComputeScore(x1e-x1s,cscore, candidateRef.getReadEndC()-candidateRef.getReadStartC(), candidateRef.getScore(),overlap);
                    cscore-=reduced_Score;

                }
                if (candidateRef.getReadStartC() < x1s) {

                    x1s = candidateRef.getReadStartC();

                }
                if (candidateRef.getReadEndC() > x1e) {

                    x1e = candidateRef.getReadEndC();

                }

            }

        }

        result.setPer_base_score(cscore/coverage);
        result.setCoveragefratction(coverageFraction);
        return result;
    }

    private static long ComputeScore(long length1, long cscore, long length2, long readscore,long overlap) {

        return (cscore/length1+readscore/length2)*overlap/2;

    }

}
