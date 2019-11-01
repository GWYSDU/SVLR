public class Parameter_Setting {

    int support_read_num=10;
    int min_sv_size=40;
    int max_sv_size=100000;
    int segment_gap_tolerance=10;
    int segment_overlap_tolerance=5;
    int partition_max_distance=5000;
    int distance_normalizer=900;
    double cluter_max_distance=0.7;
    double del_ins_dup_max_distance=1.0;
    int trans_destination_partition_max_distance=1000;
    int trans_partition_max_distance=200;
    int trans_sv_max_distance=500;
    String distance_metrics="sl";
    double min_SV_cluster_score=20;

    public Parameter_Setting(){

    }
    public Parameter_Setting(int min_sv_size,int max_sv_size,int segment_gap_tolerance, int segment_overlap_tolerance, int partition_max_distance, int distance_normalizer, double cluter_max_distance, double del_ins_dup_max_distance, int trans_destination_partition_max_distance, int trans_partition_max_distance, int trans_sv_max_distance,String distance_metrics){

        this.min_sv_size=min_sv_size;
        this.max_sv_size=max_sv_size;
        this.segment_gap_tolerance=segment_gap_tolerance;
        this.segment_overlap_tolerance=segment_overlap_tolerance;
        this.partition_max_distance=partition_max_distance;
        this.distance_normalizer=distance_normalizer;
        this.cluter_max_distance=cluter_max_distance;
        this.del_ins_dup_max_distance=del_ins_dup_max_distance;
        this.trans_destination_partition_max_distance=trans_destination_partition_max_distance;
        this.trans_partition_max_distance=trans_partition_max_distance;
        this.trans_sv_max_distance=trans_sv_max_distance;
        this.distance_metrics=distance_metrics;

    }

}
