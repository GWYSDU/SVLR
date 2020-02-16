import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

public class SV_Signature_Cluster {

    private static double calculate_score(int cluster_size, Double std_span, Double std_pos, double span, Double ave_confidence_level, int support_read_num) {

        double span_deviation_score;
        double pos_deviation_score;
        double cluster_size_weight=cluster_size/(double)support_read_num;
        double confidence_level_score;

        if(cluster_size_weight>1){

            cluster_size_weight=1;

        }

        double confidence_level_weight=ave_confidence_level/Global.per_base_best_score;
        confidence_level_score =40*confidence_level_weight+40*cluster_size_weight;

        if(cluster_size==1){

            span_deviation_score=0;
            pos_deviation_score=0;

        }
        else{

            span_deviation_score=1-Math.min(1,std_span/span);
            pos_deviation_score=1-Math.min(1,std_pos/span);

        }

        return confidence_level_score+span_deviation_score*10+pos_deviation_score*10;//

    }

    private static double calculate_score_inversion(int cluster_size, List <SignatureInversion> inv_clusters, Double std_span, Double std_pos, double span, Double ave_confidence_level, int support_read_num) {

        double span_deviation_score;
        double pos_deviation_score;

        if(cluster_size==1){

            span_deviation_score=0;
            pos_deviation_score=0;

        }
        else{

            span_deviation_score=1-Math.min(1,std_span/span);
            pos_deviation_score=1-Math.min(1,std_pos/span);

        }

        List<String> directions=new LinkedList <>();

        for(int n=0;n<inv_clusters.size();n++){

            SignatureInversion inv_cluster=inv_clusters.get(n);
            directions.add(inv_cluster.direction);

        }

        int [] directions_counts=new int[5];
        directions_counts[0]=0;
        directions_counts[1]=0;
        directions_counts[2]=0;
        directions_counts[3]=0;
        directions_counts[4]=0;

        for(int n=0;n<directions.size();n++){

            String direction=directions.get(n);

            if(direction.equals("left_fwd")){

                directions_counts[0]+=1;

            }
            if(direction.equals("left_rev")){

                directions_counts[1]+=1;

            }
            if(direction.equals("right_fwd")){

                directions_counts[2]+=1;

            }
            if(direction.equals("right_rev")){

                directions_counts[3]+=1;

            }
            if(direction.equals("all")){

                directions_counts[4]+=1;

            }

        }

        int left_signatures=directions_counts[0]+directions_counts[1];
        int ritght_signatures=directions_counts[2]+directions_counts[3];
        int valid_suppl_signatures= Math.min(left_signatures,ritght_signatures)+directions_counts[4];
        double cluster_size_weight=valid_suppl_signatures/(double)support_read_num;

        if(cluster_size_weight>1){

            cluster_size_weight=1;

        }
        double confidence_level_score;
        double confidence_level_weight=ave_confidence_level/Global.per_base_best_score;
        confidence_level_score=40*confidence_level_weight+40*cluster_size_weight;
        return confidence_level_score+span_deviation_score*10+pos_deviation_score*10;

    }

    public static List cluster_sv_signature(Parameter_Setting Parameter_Setting) throws FileNotFoundException {

        List clusters=new LinkedList();
        clusters.add(deletion_partition_and_cluster_unilocal(Global.deletion_signatures, Parameter_Setting));
        clusters.add(insertion_partition_and_cluster_unilocal(Global.insertion_signatures, Parameter_Setting));
        clusters.add(inversion_partition_and_cluster_unilocal(Global.inversion_signatures, Parameter_Setting));
        clusters.add(tandem_duplication_partition_and_cluster_bilocal(Global.duplication_tandem_signatures, Parameter_Setting));
        clusters.add(potential_interspersed_trans_partition_and_cluster_bilocal(Global.Potential_Trans_Inters_signatures, Parameter_Setting));
        return clusters;
    }

    private static List<SignatureClusterUniLocal> deletion_partition_and_cluster_unilocal(List <SignatureDeletion> deletion_signatures, Parameter_Setting Parameter_Setting) throws FileNotFoundException {

        List<SignatureClusterUniLocal> sort_deletion_consolidate_clusters_unilocal=new LinkedList <>();

        if(deletion_signatures!=null){

            if(deletion_signatures.size()>0){

                List<List<SignatureDeletion>> partitions=deletion_form_partitions(deletion_signatures,Parameter_Setting.partition_max_distance);
                List<List<SignatureDeletion>> clusters=deletion_clusters_from_partitions(partitions,Parameter_Setting);
                sort_deletion_consolidate_clusters_unilocal=deletion_consolidate_clusters_unilocal(clusters,Parameter_Setting);
                Collections.sort(sort_deletion_consolidate_clusters_unilocal,new contig_position());

            }

        }

        return sort_deletion_consolidate_clusters_unilocal;

    }

    private static List <List <SignatureDeletion>> deletion_form_partitions(List <SignatureDeletion> signatures, int partition_max_distance) {

        List<List<SignatureDeletion>> partitions=new LinkedList<>();
        List<SignatureDeletion> current_partition=new LinkedList <>();
        int signaturesize=signatures.size();
        Collections.sort(signatures,new deletion_key_sort());

        for(int n=0;n<signaturesize;n++){

            SignatureDeletion signatureDeletion=signatures.get(n);

            if(current_partition.size()>0&&current_partition.get(0).mean_distance_to(signatureDeletion)>partition_max_distance){

                partitions.add(current_partition);
                current_partition=new LinkedList <>();

            }

            current_partition.add(signatureDeletion);

        }

        if(current_partition.size()>0){

            partitions.add(current_partition);

        }

        return partitions;

    }

    private static class deletion_key_sort implements Comparator<SignatureDeletion> {


        @Override
        public int compare(SignatureDeletion o1, SignatureDeletion o2) {

           if(o1.get_key_contig().length()>o2.get_key_contig().length()){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())>0){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()>o2.get_key_position()){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()==o2.get_key_position()){

                return 0;

            }
            else{

                return -1;

            }

        }

    }

    private static List <List <SignatureDeletion>> deletion_clusters_from_partitions(List<List<SignatureDeletion>> partitions, Parameter_Setting Parameter_Setting) {

        List<List<SignatureDeletion>> clusters_full=new LinkedList <>();
        int partition_size=partitions.size();

        for(int n=0;n<partition_size;n++){

            List<SignatureDeletion> partition_sample=new ArrayList <>();
            List<SignatureDeletion> partition=partitions.get(n);
            int i;
            Random random=new Random();

            if(partition.size()>100){

                while (true){

                    if(partition_sample.size()>100){

                        break;

                    }
                    i=random.nextInt(partition.size());
                    partition_sample.add(partition.get(i));

                }

            }
            else{

                partition_sample=partition;

            }

            Collections.sort(partition_sample,new deletion_indel_size());
            Global.cliques=new ArrayList <>();
            MaximalCliquesWithPivot connection_graph=new MaximalCliquesWithPivot(partition_sample.size());
            connection_graph.Deletion_creatGraph(partition_sample,Parameter_Setting);
            connection_graph.Bron_KerboschpivotExecute();
            int cliques_size=Global.cliques.size();

            if(cliques_size>0){

                for(int j=0;j<cliques_size;j++){

                    List<Vertex> clique=Global.cliques.get(j);

                    int clique_size=clique.size();

                    List<SignatureDeletion> cluster=new LinkedList <>();

                    for(int m=0;m<clique_size;m++){

                        cluster.add(partition_sample.get(clique.get(m).getX()));

                    }

                    clusters_full.add(cluster);

                }

            }

        }

        return clusters_full;

    }

    private static class deletion_indel_size implements Comparator <SignatureDeletion> {

        @Override
        public int compare(SignatureDeletion o1, SignatureDeletion o2) {

            long indel_size1=o1.getSourceEnd()-o1.getSourceStart();
            long indel_size2=o2.getSourceEnd()-o2.getSourceStart();

            if(indel_size1>indel_size2){

                return  1;

            }
            else if(indel_size1==indel_size2){

                return  0;

            }
            else{

                return -1;

            }

        }

    }

    private static List <SignatureClusterUniLocal> deletion_consolidate_clusters_unilocal(List<List<SignatureDeletion>> clusters, Parameter_Setting options) {

        List<SignatureClusterUniLocal> consolidated_cluters=new LinkedList <>();

        if(clusters!=null){

            int clusters_size=clusters.size();

            if(clusters_size>0){

                for(int n=0;n<clusters_size;n++){

                    List<SignatureDeletion> cluster=clusters.get(n);
                    int cluster_size=cluster.size();

                    if(cluster_size>0){

                        SignatureDeletion signatureDeletion=cluster.get(0);
                        Double std_span;
                        Double std_pos;
                        double score;
                        long start_sum=0;
                        long end_sum=0;
                        double ave_span=0;
                        double confidence_level_sum=0;
                        SignatureDeletion member;

                        for(int m=0;m<cluster_size;m++){

                            member=cluster.get(m);
                            start_sum+=member.getSourceStart();
                            end_sum+=member.getSourceEnd();
                            confidence_level_sum+=member.getConfidence_level();

                        }

                        double average_start=(double) start_sum/cluster_size;
                        double average_end=(double)end_sum/cluster_size;
                        double average_confidence_level=confidence_level_sum/cluster_size;

                        if(cluster_size>1){

                            long span=0;
                            long pos=0;

                            for(int m=0;m<cluster_size;m++){

                                member=cluster.get(m);
                                span+=member.span;
                                pos+=(member.getSourceEnd()+member.getSourceStart())/2;

                            }

                            ave_span=(double)span/cluster_size;
                            double ave_pos=(double)pos/cluster_size;

                            double span_qua_sum=0;
                            double pos_qua_sum=0;

                            for(int m=0;m<cluster_size;m++){

                                member=cluster.get(m);
                                span_qua_sum+=Math.pow(member.span-ave_span,2);
                                pos_qua_sum+=Math.pow((member.getSourceEnd()+member.getSourceStart())/2-ave_pos,2);

                            }

                            std_span=Math.sqrt(span_qua_sum/cluster_size);
                            std_pos=Math.sqrt(pos_qua_sum/cluster_size);

                        }
                        else{

                            std_span=0.0;
                            std_pos= 0.0;

                        }

                        score=calculate_score(cluster_size,std_span,std_pos,average_end-average_start,average_confidence_level,options.support_read_num);

                        if(score>=options.min_SV_cluster_score&&(average_end-average_start)<150000){

                            consolidated_cluters.add(new SignatureClusterUniLocal(signatureDeletion.getSourceContig(),(int)(Math.round(average_start)),(int)(Math.round(average_end)),score,cluster.size(),cluster,signatureDeletion.getType(),std_span,std_pos,ave_span));

                        }

                    }

                }

            }

        }

        return consolidated_cluters;

    }

    private static List<SignatureClusterUniLocal>  insertion_partition_and_cluster_unilocal(List <SignatureInsertion> insertion_signatures, Parameter_Setting options) {

        List<SignatureClusterUniLocal> sort_insertion_consolidate_clusters_unilocal=new LinkedList <>();
        if(insertion_signatures!=null){

            if(insertion_signatures.size()>0){

                List <List <SignatureInsertion>> partitions = insertion_form_partitions(insertion_signatures, options.partition_max_distance);
                List <List <SignatureInsertion>> clusters = insertion_clusters_from_partitions(partitions, options);
                sort_insertion_consolidate_clusters_unilocal=insertion_consolidate_clusters_unilocal(clusters, options);
                Collections.sort(sort_insertion_consolidate_clusters_unilocal,new contig_position());

            }

        }

        return sort_insertion_consolidate_clusters_unilocal;
    }

    private static List<List<SignatureInsertion>> insertion_form_partitions(List<SignatureInsertion> insertion_signatures, int partition_max_distance) {

        List<List<SignatureInsertion>> partitions=new LinkedList <>();
        List<SignatureInsertion> current_partition=new LinkedList <>();
        int signaturesize=insertion_signatures.size();
        Collections.sort(insertion_signatures,new insertion_key_sort());

        for(int n=0;n<signaturesize;n++){

            SignatureInsertion SignatureInsertion=insertion_signatures.get(n);

            if(current_partition.size()>0&&current_partition.get(0).mean_distance_to(SignatureInsertion)>partition_max_distance){

                partitions.add(current_partition);
                current_partition=new LinkedList <>();

            }

            current_partition.add(SignatureInsertion);

        }

        if(current_partition.size()>0){

            partitions.add(current_partition);

        }

        return partitions;

    }

    private static class insertion_key_sort implements Comparator<SignatureInsertion>{

        @Override
        public int compare(SignatureInsertion o1, SignatureInsertion o2) {

            if(o1.get_key_contig().length()>o2.get_key_contig().length()){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())>0){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()>o2.get_key_position()){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()==o2.get_key_position()){

                return 0;

            }
            else{

                return -1;

            }

        }

    }

    private static List<List<SignatureInsertion>> insertion_clusters_from_partitions(List<List<SignatureInsertion>> partitions, Parameter_Setting options) {

        List<List<SignatureInsertion>> clusters_full=new LinkedList <>();
        int partition_size=partitions.size();

        for(int n=0;n<partition_size;n++){

            List<SignatureInsertion> partition_sample=new ArrayList <>();
            List<SignatureInsertion> partition=partitions.get(n);
            int i;
            Random random=new Random();

            if(partition.size()>100){

                while (true){

                    if(partition_sample.size()>100){

                        break;

                    }
                    i=random.nextInt(partition.size());
                    partition_sample.add(partition.get(i));

                }

            }
            else{

                partition_sample=partition;

            }

            Collections.sort(partition_sample,new insertion_indel_size());
            SignatureInsertion largest_signature=partition_sample.get(partition_sample.size()-1);
            long largest_indel_size=largest_signature.getSourceEnd()-largest_signature.getSourceStart();
            Global.cliques=new ArrayList <>();
            MaximalCliquesWithPivot connection_graph=new MaximalCliquesWithPivot(partition_sample.size());
            connection_graph.Insertion_creatGraph(partition_sample,options,largest_indel_size);
            connection_graph.Bron_KerboschpivotExecute();
            int cliques_size=Global.cliques.size();

            if(cliques_size>0){

                for(int j=0;j<cliques_size;j++){

                    List<Vertex> clique=Global.cliques.get(j);

                    int clique_size=clique.size();

                    List<SignatureInsertion> cluster=new LinkedList <>();

                    for(int m=0;m<clique_size;m++){

                        cluster.add(partition_sample.get(clique.get(m).getX()));

                    }

                    clusters_full.add(cluster);

                }

            }

        }

        return clusters_full;

    }

    private static class insertion_indel_size implements Comparator<SignatureInsertion>{

        @Override
        public int compare(SignatureInsertion o1, SignatureInsertion o2) {

            long indel_size1=o1.getSourceEnd()-o1.getSourceStart();
            long indel_size2=o2.getSourceEnd()-o2.getSourceStart();

            if(indel_size1>indel_size2){

                return  1;

            }
            else if(indel_size1==indel_size2){

                return  0;

            }
            else{

                return -1;

            }

        }

    }

    private static List<SignatureClusterUniLocal> insertion_consolidate_clusters_unilocal(List<List<SignatureInsertion>> clusters, Parameter_Setting options) {

        List<SignatureClusterUniLocal> consolidated_cluters=new LinkedList <>();

        if(clusters!=null){

            int clusters_size=clusters.size();

            for(int n=0;n<clusters_size;n++){

                List<SignatureInsertion> cluster=clusters.get(n);
                int cluster_size=cluster.size();

                if(cluster_size>0){

                    SignatureInsertion SignatureInsertion=cluster.get(0);
                    Double std_span;
                    Double std_pos;
                    double score;
                    long start_sum=0;
                    long end_sum=0;
                    double ave_span=0;
                    double confidence_level_sum=0;
                    SignatureInsertion member;

                    for(int m=0;m<cluster_size;m++){

                        member=cluster.get(m);
                        start_sum+=member.getSourceStart();
                        end_sum+=member.getSourceEnd();
                        confidence_level_sum+=member.getConfidence_level();

                    }

                    double average_start=(double) start_sum/cluster_size;
                    double average_end=(double)end_sum/cluster_size;
                    double average_confidence_level=confidence_level_sum/cluster_size;

                    if(cluster_size>1){

                        long span=0;
                        long pos=0;

                        for(int m=0;m<cluster_size;m++){

                            member=cluster.get(m);
                            span+=member.span;
                            pos+=(member.getSourceEnd()+member.getSourceStart())/2;

                        }

                        ave_span=(double)span/cluster_size;
                        double ave_pos=(double)pos/cluster_size;

                        double span_qua_sum=0;
                        double pos_qua_sum=0;

                        for(int m=0;m<cluster_size;m++){

                            member=cluster.get(m);
                            span_qua_sum+=Math.pow(member.getSourceEnd()-member.getSourceStart()-ave_span,2);
                            pos_qua_sum+=Math.pow((member.getSourceEnd()+member.getSourceStart())/2-ave_pos,2);

                        }

                        std_span=Math.sqrt(span_qua_sum/cluster_size);
                        std_pos=Math.sqrt(pos_qua_sum/cluster_size);

                    }
                    else{

                        std_span= 0.0;
                        std_pos= 0.0;

                    }

                    score=calculate_score(cluster_size, std_span,std_pos,average_end-average_start,average_confidence_level, options.support_read_num);

                    if(score>=options.min_SV_cluster_score){

                        consolidated_cluters.add(new SignatureClusterUniLocal(SignatureInsertion.getSourceContig(),(int)(Math.round(average_start)),(int)(Math.round(average_end)),score,cluster.size(),cluster,SignatureInsertion.getType(),std_span,std_pos,ave_span));

                    }

                }

            }

        }

        return consolidated_cluters;

    }

    private static List<SignatureClusterUniLocal> inversion_partition_and_cluster_unilocal(List <SignatureInversion> inversion_signatures, Parameter_Setting options) {

        List<SignatureClusterUniLocal> sort_inversion_consolidate_clusters_unilocal=new LinkedList <>();
        if(inversion_signatures!=null){

            if(inversion_signatures.size()>0){

                List<List<SignatureInversion>> partitions=inversion_form_partitions(inversion_signatures,options.partition_max_distance);
                List<List<SignatureInversion>> clusters=inversion_clusters_from_partitions(partitions,options);
                sort_inversion_consolidate_clusters_unilocal =inversion_consolidate_clusters_unilocal(clusters,options);
                Collections.sort(sort_inversion_consolidate_clusters_unilocal,new contig_position());

            }

        }

        return sort_inversion_consolidate_clusters_unilocal;

    }

    private static List <List <SignatureInversion>> inversion_form_partitions(List <SignatureInversion> signatures, int partition_max_distance) {

        List<List<SignatureInversion>> partitions=new LinkedList <>();
        List<SignatureInversion> current_partition=new LinkedList <>();
        int signaturesize=signatures.size();
        Collections.sort(signatures,new inversion_key_sort());

        for(int n=0;n<signaturesize;n++){

            SignatureInversion SignatureInversion=signatures.get(n);

            if(current_partition.size()>0&&current_partition.get(0).mean_distance_to(SignatureInversion)>partition_max_distance){

                partitions.add(current_partition);
                current_partition=new LinkedList <>();

            }

            current_partition.add(SignatureInversion);

        }

        if(current_partition.size()>0){

            partitions.add(current_partition);

        }

        return partitions;

    }

    private static class inversion_key_sort implements Comparator <SignatureInversion> {


        @Override
        public int compare(SignatureInversion o1, SignatureInversion o2) {

            if(o1.get_key_contig().length()>o2.get_key_contig().length()){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())>0){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()>o2.get_key_position()){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()==o2.get_key_position()){

                return 0;

            }
            else{

                return -1;

            }

        }

    }

    private static List <List <SignatureInversion>> inversion_clusters_from_partitions(List<List<SignatureInversion>> partitions, Parameter_Setting options) {

        List<List<SignatureInversion>> clusters_full=new LinkedList <>();
        int partition_size=partitions.size();

        for(int n=0;n<partition_size;n++){

            List<SignatureInversion> partition_sample=new ArrayList <>();
            List<SignatureInversion> partition=partitions.get(n);
            int i;
            Random random=new Random();

            if(partition.size()>100){

                while (true){

                    if(partition_sample.size()>100){

                        break;

                    }
                    i=random.nextInt(partition.size());
                    partition_sample.add(partition.get(i));

                }

            }
            else{

                partition_sample= partition;

            }

            Collections.sort(partition_sample,new inversion_indel_size());
            SignatureInversion largest_signature=partition_sample.get(partition_sample.size()-1);
            long largest_indel_size=largest_signature.getSourceEnd()-largest_signature.getSourceStart();
            Global.cliques=new ArrayList <>();
            MaximalCliquesWithPivot connection_graph=new MaximalCliquesWithPivot(partition_sample.size());
            connection_graph.Inversion_creatGraph(partition_sample,options,largest_indel_size);
            connection_graph.Bron_KerboschpivotExecute();
            int cliques_size=Global.cliques.size();

            if(cliques_size>0){

                for(int j=0;j<cliques_size;j++){

                    List<Vertex> clique=Global.cliques.get(j);

                    int clique_size=clique.size();

                    List<SignatureInversion> cluster=new LinkedList <>();

                    for(int m=0;m<clique_size;m++){

                        cluster.add(partition_sample.get(clique.get(m).getX()));

                    }

                    clusters_full.add(cluster);

                }

            }

        }

        return clusters_full;

    }

    private static class inversion_indel_size implements Comparator <SignatureInversion> {

        @Override
        public int compare(SignatureInversion o1, SignatureInversion o2) {

            long indel_size1=o1.getSourceEnd()-o1.getSourceStart();
            long indel_size2=o2.getSourceEnd()-o2.getSourceStart();

            if(indel_size1>indel_size2){

                return  1;

            }
            else if(indel_size1==indel_size2){

                return  0;

            }
            else{

                return -1;

            }

        }

    }


    private static List <SignatureClusterUniLocal> inversion_consolidate_clusters_unilocal(List<List<SignatureInversion>> clusters, Parameter_Setting options) {

        List<SignatureClusterUniLocal> consolidated_cluters=new LinkedList <>();

        if(clusters!=null){

            int clusters_size=clusters.size();

            for(int n=0;n<clusters_size;n++){

                List<SignatureInversion> cluster=clusters.get(n);
                int cluster_size=cluster.size();

                if(cluster_size>0){

                    SignatureInversion SignatureInversion=cluster.get(0);
                    Double std_span;
                    Double std_pos;
                    double score;
                    long start_sum=0;
                    long end_sum=0;
                    double ave_span=0;
                    double confidence_level_sum=0;
                    SignatureInversion member;

                    for(int m=0;m<cluster_size;m++){

                        member=cluster.get(m);
                        start_sum+=member.getSourceStart();
                        end_sum+=member.getSourceEnd();
                        confidence_level_sum+=member.getConfidence_level();

                    }

                    double average_start=(double) start_sum/cluster_size;
                    double average_end=(double)end_sum/cluster_size;
                    double average_confidence_level=confidence_level_sum/cluster_size;

                    if(cluster_size>1){

                        long span=0;
                        long pos=0;

                        for(int m=0;m<cluster_size;m++){

                            member=cluster.get(m);
                            span+=member.getSourceEnd()-member.getSourceStart();
                            pos+=(member.getSourceEnd()+member.getSourceStart())/2;

                        }

                        ave_span=(double)span/cluster_size;
                        double ave_pos=(double)pos/cluster_size;

                        double span_qua_sum=0;
                        double pos_qua_sum=0;

                        for(int m=0;m<cluster_size;m++){

                            member=cluster.get(m);
                            span_qua_sum+=Math.pow(member.getSourceEnd()-member.getSourceStart()-ave_span,2);
                            pos_qua_sum+=Math.pow((member.getSourceEnd()+member.getSourceStart())/2-ave_pos,2);

                        }

                        std_span=Math.sqrt(span_qua_sum/cluster_size);
                        std_pos=Math.sqrt(pos_qua_sum/cluster_size);

                    }
                    else{

                        std_span= 0.0;
                        std_pos= 0.0;

                    }

                    score=calculate_score_inversion(cluster_size,cluster,std_span,std_pos,average_end-average_start,average_confidence_level,options.support_read_num);

                    if(score>=options.min_SV_cluster_score){

                        consolidated_cluters.add(new SignatureClusterUniLocal(SignatureInversion.getSourceContig(),(int)(Math.round(average_start)),(int)(Math.round(average_end)),score,cluster.size(),cluster,SignatureInversion.getType(),std_span,std_pos,ave_span));

                    }

                }

            }

        }

        return consolidated_cluters;

    }

    private static List<SignatureClusterBilocal> tandem_duplication_partition_and_cluster_bilocal(List <SignatureDuplicationTandem> tandem_duplication_signatures, Parameter_Setting options) {

        List<SignatureClusterBilocal> sort_duplication_tandem_consolidate_clusters_bilocal=new LinkedList <>();
        if(tandem_duplication_signatures!=null){

            if(tandem_duplication_signatures.size()>0){

                List<List<SignatureDuplicationTandem>> partitions=duplication_tandem_form_partitions(tandem_duplication_signatures,options.partition_max_distance);
                List<List<SignatureDuplicationTandem>> clusters=duplication_tandem_clusters_from_partitions(partitions,options);
                sort_duplication_tandem_consolidate_clusters_bilocal=duplication_tandem_consolidate_clusters_bilocal(clusters,options);
                Collections.sort(sort_duplication_tandem_consolidate_clusters_bilocal,new bilocal_contig_position());
            }

        }

        return sort_duplication_tandem_consolidate_clusters_bilocal;

    }

    private static List <List <SignatureDuplicationTandem>> duplication_tandem_form_partitions(List <SignatureDuplicationTandem> signatures, int partition_max_distance) {

        List<List<SignatureDuplicationTandem>> partitions=new LinkedList <>();
        List<SignatureDuplicationTandem> current_partition=new LinkedList <>();
        int signaturesize=signatures.size();
        Collections.sort(signatures,new duplication_tandem_key_sort());

        for(int n=0;n<signaturesize;n++){

            SignatureDuplicationTandem SignatureDuplicationTandem=signatures.get(n);

            if(current_partition.size()>0&&current_partition.get(0).mean_distance_to(SignatureDuplicationTandem)>partition_max_distance){

                partitions.add(current_partition);
                current_partition=new LinkedList <>();

            }

            current_partition.add(SignatureDuplicationTandem);

        }

        if(current_partition.size()>0){

            partitions.add(current_partition);

        }

        return partitions;

    }

    private static class duplication_tandem_key_sort implements Comparator <SignatureDuplicationTandem> {


        @Override
        public int compare(SignatureDuplicationTandem o1, SignatureDuplicationTandem o2) {

            if(o1.get_key_contig().length()>o2.get_key_contig().length()){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())>0){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()>o2.get_key_position()){

                return 1;

            }
            else if(o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()==o2.get_key_position()){

                return 0;

            }
            else{

                return -1;

            }

        }

    }

    private static List <List <SignatureDuplicationTandem>> duplication_tandem_clusters_from_partitions(List<List<SignatureDuplicationTandem>> partitions, Parameter_Setting options) {

        List<List<SignatureDuplicationTandem>> clusters_full=new LinkedList <>();
        int partition_size=partitions.size();

        for(int n=0;n<partition_size;n++){

            List<SignatureDuplicationTandem> partition_sample=new ArrayList <>();
            List<SignatureDuplicationTandem> partition=partitions.get(n);
            int i;
            Random random=new Random();

            if(partition.size()>100){

                while (true){

                    if(partition_sample.size()>100){

                        break;

                    }
                    i=random.nextInt(partition.size());
                    partition_sample.add(partition.get(i));

                }

            }
            else{

                partition_sample=partition;

            }

            Collections.sort(partition_sample,new duplication_tandem_indel_size());
            Global.cliques=new ArrayList <>();
            MaximalCliquesWithPivot connection_graph=new MaximalCliquesWithPivot(partition_sample.size());
            connection_graph.Duplication_Tandem_creatGraph(partition_sample,options);
            connection_graph.Bron_KerboschpivotExecute();
            int cliques_size=Global.cliques.size();

            if(cliques_size>0){

                for(int j=0;j<cliques_size;j++){

                    List<Vertex> clique=Global.cliques.get(j);

                    int clique_size=clique.size();

                    List<SignatureDuplicationTandem> cluster=new LinkedList <>();

                    for(int m=0;m<clique_size;m++){

                        cluster.add(partition_sample.get(clique.get(m).getX()));

                    }

                    clusters_full.add(cluster);

                }

            }

        }

        return clusters_full;

    }

    private static class duplication_tandem_indel_size implements Comparator <SignatureDuplicationTandem> {

        @Override
        public int compare(SignatureDuplicationTandem o1, SignatureDuplicationTandem o2) {

            long indel_size1=o1.getSourceEnd()-o1.getSourceStart();
            long indel_size2=o2.getSourceEnd()-o2.getSourceStart();

            if(indel_size1>indel_size2){

                return  1;

            }
            else if(indel_size1==indel_size2){

                return  0;

            }
            else{

                return -1;

            }

        }

    }


    private static List <SignatureClusterBilocal> duplication_tandem_consolidate_clusters_bilocal(List<List<SignatureDuplicationTandem>> clusters, Parameter_Setting options) {

        List<SignatureClusterBilocal> consolidated_cluters=new LinkedList <>();

        if(clusters!=null){

            int clusters_size=clusters.size();

            for(int n=0;n<clusters_size;n++){

                List<SignatureDuplicationTandem> cluster=clusters.get(n);
                int cluster_size=cluster.size();

                if(cluster_size>0){

                    SignatureDuplicationTandem SignatureDuplicationTandem=cluster.get(0);
                    long source_start_sum=0;
                    long source_end_sum=0;
                    double confidence_level_sum=0;
                    SignatureDuplicationTandem member;

                    Double source_std_span;
                    Double source_std_pos;
                    double score;

                    for(int m=0;m<cluster_size;m++){

                        member=cluster.get(m);
                        source_start_sum+=member.getSourceStart();
                        source_end_sum+=member.getSourceEnd();
                        confidence_level_sum+=member.getConfidence_level();

                    }

                    double source_average_start=(double) source_start_sum/cluster_size;
                    double source_average_end=(double)source_end_sum/cluster_size;
                    double average_confidence_level=confidence_level_sum/cluster_size;

                    if(cluster_size>1){

                        long source_span=0;
                        long source_pos=0;

                        for(int m=0;m<cluster_size;m++){

                            member=cluster.get(m);
                            source_span+=member.getSourceEnd()-member.getSourceStart();
                            source_pos+=(member.getSourceEnd()+member.getSourceStart())/2;

                        }

                        double source_ave_span=(double)source_span/cluster_size;
                        double source_ave_pos=(double)source_pos/cluster_size;

                        double source_span_qua_sum=0;
                        double source_pos_qua_sum=0;

                        for(int m=0;m<cluster_size;m++){

                            member=cluster.get(m);
                            source_span_qua_sum+=Math.pow(member.getSourceEnd()-member.getSourceStart()-source_ave_span,2);
                            source_pos_qua_sum+=Math.pow((member.getSourceEnd()+member.getSourceStart())/2-source_ave_pos,2);

                        }

                        source_std_span=Math.sqrt(source_span_qua_sum/cluster_size);
                        source_std_pos=Math.sqrt(source_pos_qua_sum/cluster_size);

                    }
                    else{

                        source_std_span= 0.0;
                        source_std_pos=0.0;

                    }

                    int max_copies=0;

                    for(int i=0;cluster.size() > i; i++){

                        int copies=cluster.get(i).copies;

                        if(max_copies<=copies){

                            max_copies=copies;

                        }

                    }

                    score=calculate_score(cluster_size, source_std_span,source_std_pos,source_average_end-source_average_start,average_confidence_level, options.support_read_num);

                    if(score>options.min_SV_cluster_score&&max_copies>1){

                        consolidated_cluters.add(new SignatureClusterBilocal(SignatureDuplicationTandem.getSourceContig(),(int)(Math.round(source_average_start)),(int)(Math.round(source_average_end)),SignatureDuplicationTandem.getSourceContig(),(int)(Math.round(source_average_end)),(int)(Math.round(source_average_end))+max_copies*((int)(Math.round(source_average_end))-(int)(Math.round(source_average_start))),score,cluster.size(),cluster,SignatureDuplicationTandem.getType(),source_std_span,source_std_pos));

                    }

                }

            }

        }

        return consolidated_cluters;

    }

    private static List<SignatureClusterBilocal> potential_interspersed_trans_partition_and_cluster_bilocal(List <SignaturePotential_Trans_Inters> potential_interspersed_trans_signatures, Parameter_Setting options) {

        List<SignatureClusterBilocal> sort_interspersedduplication_consolidate_clusters_bilocal=new LinkedList <>();

        if(potential_interspersed_trans_signatures!=null){

            if(potential_interspersed_trans_signatures.size()>0){

                List<List<SignaturePotential_Trans_Inters>> partitions= potential_interspersed_trans_form_partitions(potential_interspersed_trans_signatures,options.partition_max_distance);
                List<List<SignaturePotential_Trans_Inters>> clusters= potential_interspersed_trans_clusters_from_partitions(partitions,options);
                sort_interspersedduplication_consolidate_clusters_bilocal= potential_interspersed_trans_consolidate_clusters_bilocal(clusters,options);
                Collections.sort(sort_interspersedduplication_consolidate_clusters_bilocal,new bilocal_contig_position());

            }

        }

        return sort_interspersedduplication_consolidate_clusters_bilocal;

    }

    private static List <List <SignaturePotential_Trans_Inters>> potential_interspersed_trans_form_partitions(List <SignaturePotential_Trans_Inters> signatures, int partition_max_distance) {

        List<List<SignaturePotential_Trans_Inters>> partitions=new LinkedList <>();
        List<SignaturePotential_Trans_Inters> current_partition=new LinkedList <>();
        int signaturesize=signatures.size();
        Collections.sort(signatures,new potential_interspersed_trans_key_sort());

        for(int n=0;n<signaturesize;n++){

            SignaturePotential_Trans_Inters SignaturePotential_Trans_Inters=signatures.get(n);

            if(current_partition.size()>0&&current_partition.get(0).mean_distance_to(SignaturePotential_Trans_Inters)>partition_max_distance){

                partitions.add(current_partition);
                current_partition=new LinkedList <>();

            }

            current_partition.add(SignaturePotential_Trans_Inters);

        }

        if(current_partition.size()>0){

            partitions.add(current_partition);

        }

        return partitions;

    }

    private static class potential_interspersed_trans_key_sort implements Comparator <SignaturePotential_Trans_Inters> {


        @Override
        public int compare(SignaturePotential_Trans_Inters o1, SignaturePotential_Trans_Inters o2) {

           if(o1.get_key_contig().length()>o2.get_key_contig().length()){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())>0){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()>o2.get_key_position()){

                return 1;

            }
            else if(o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()==o2.get_key_position()){

                return 0;

            }
            else{

                return -1;

            }

        }

    }

    private static List <List <SignaturePotential_Trans_Inters>> potential_interspersed_trans_clusters_from_partitions(List <List <SignaturePotential_Trans_Inters>> partitions, Parameter_Setting options) {

        List<List<SignaturePotential_Trans_Inters>> clusters_full=new LinkedList <>();
        int partition_size=partitions.size();

        for(int n=0;n<partition_size;n++){

            List<SignaturePotential_Trans_Inters> partition_sample=new ArrayList <>();
            List<SignaturePotential_Trans_Inters> partition=partitions.get(n);
            int i;
            Random random=new Random();

            if(partition.size()>100){

                while (true){

                    if(partition_sample.size()>100){

                        break;

                    }
                    i=random.nextInt(partition.size());
                    partition_sample.add(partition.get(i));

                }

            }
            else{

                partition_sample= partition;

            }

            Collections.sort(partition_sample,new potential_interspersed_trans_indel_size());
            Global.cliques=new ArrayList <>();
            MaximalCliquesWithPivot connection_graph=new MaximalCliquesWithPivot(partition_sample.size());
            connection_graph.potential_interspersed_trans_creatGraph(partition_sample,options);
            connection_graph.Bron_KerboschpivotExecute();
            int cliques_size=Global.cliques.size();

            if(cliques_size>0){

                for(int j=0;j<cliques_size;j++){

                    List<Vertex> clique=Global.cliques.get(j);

                    int clique_size=clique.size();

                    List<SignaturePotential_Trans_Inters> cluster=new LinkedList <>();

                    for(int m=0;m<clique_size;m++){

                        cluster.add(partition_sample.get(clique.get(m).getX()));

                    }

                    clusters_full.add(cluster);

                }

            }

        }

        return clusters_full;

    }

    private static class potential_interspersed_trans_indel_size implements Comparator <SignaturePotential_Trans_Inters> {

        @Override
        public int compare(SignaturePotential_Trans_Inters o1, SignaturePotential_Trans_Inters o2) {

            long indel_size1=o1.getSourceEnd()-o1.getSourceStart();
            long indel_size2=o2.getSourceEnd()-o2.getSourceStart();

            if(indel_size1>indel_size2){

                return  1;

            }
            else if(indel_size1==indel_size2){

                return  0;

            }
            else{

                return -1;

            }

        }

    }


    private static List <SignatureClusterBilocal> potential_interspersed_trans_consolidate_clusters_bilocal(List<List<SignaturePotential_Trans_Inters>> clusters, Parameter_Setting options) {

        List<SignatureClusterBilocal> consolidated_cluters=new LinkedList <>();

        if(clusters!=null){

            int clusters_size=clusters.size();

            for(int n=0;n<clusters_size;n++){

                List<SignaturePotential_Trans_Inters> cluster=clusters.get(n);
                int cluster_size=cluster.size();

                if(cluster_size>0){

                    SignaturePotential_Trans_Inters SignaturePotential_Trans_Inters=cluster.get(0);
                    long source_start_sum=0;
                    long source_end_sum=0;
                    double confidence_level_sum=0;
                    SignaturePotential_Trans_Inters member;

                    Double source_std_span;
                    Double source_std_pos;
                    double score;

                    for(int m=0;m<cluster_size;m++){

                        member=cluster.get(m);
                        source_start_sum+=member.getSourceStart();
                        source_end_sum+=member.getSourceEnd();
                        confidence_level_sum+=member.getConfidence_level();

                    }

                    double source_average_start=(double) source_start_sum/cluster_size;
                    double source_average_end=(double)source_end_sum/cluster_size;
                    double average_confidence_level=confidence_level_sum/cluster_size;

                    if(cluster_size>1){

                        long source_span=0;
                        long source_pos=0;

                        for(int m=0;m<cluster_size;m++){

                            member=cluster.get(m);
                            source_span+=member.getSourceEnd()-member.getSourceStart();
                            source_pos+=(member.getSourceEnd()+member.getSourceStart())/2;

                        }

                        double source_ave_span=(double)source_span/cluster_size;
                        double source_ave_pos=(double)source_pos/cluster_size;

                        double source_span_qua_sum=0;
                        double source_pos_qua_sum=0;

                        for(int m=0;m<cluster_size;m++){

                            member=cluster.get(m);
                            source_span_qua_sum+=Math.pow(member.getSourceEnd()-member.getSourceStart()-source_ave_span,2);
                            source_pos_qua_sum+=Math.pow((member.getSourceEnd()+member.getSourceStart())/2-source_ave_pos,2);

                        }

                        source_std_span=Math.sqrt(source_span_qua_sum/cluster_size);
                        source_std_pos=Math.sqrt(source_pos_qua_sum/cluster_size);

                    }
                    else{

                        source_std_span= 0.0;
                        source_std_pos= 0.0;

                    }

                    long destination_start_sum=0;
                    long destination_end_sum=0;
                    Double destination_std_span;
                    Double destination_std_pos;

                    for(int m=0;m<cluster_size;m++){

                        destination_start_sum+=cluster.get(m).getDestinationStart();

                    }

                    double destination_average_start=(double) destination_start_sum/cluster_size;

                    for(int m=0;m<cluster_size;m++){

                        destination_end_sum+=cluster.get(m).getDestinationEnd();

                    }

                    double destination_average_end=(double)destination_end_sum/cluster_size;

                    if(cluster_size>1){

                        long destination_span=0;
                        long destination_pos=0;

                        for(int m=0;m<cluster_size;m++){

                            member=cluster.get(m);
                            destination_span+=member.getDestinationEnd()-member.getDestinationStart();
                            destination_pos+=(member.getDestinationEnd()+member.getDestinationStart())/2;

                        }

                        double destination_ave_span=(double)destination_span/cluster_size;
                        double destination_ave_pos=(double)destination_pos/cluster_size;

                        double destination_span_qua_sum=0;
                        double destination_pos_qua_sum=0;

                        for(int m=0;m<cluster_size;m++){

                            member=cluster.get(m);
                            destination_span_qua_sum+=Math.pow(member.getDestinationEnd()-member.getDestinationStart()-destination_ave_span,2);
                            destination_pos_qua_sum+=Math.pow((member.getDestinationEnd()+member.getDestinationStart())/2-destination_ave_pos,2);

                        }

                        destination_std_span=Math.sqrt(destination_span_qua_sum/cluster_size);
                        destination_std_pos=Math.sqrt(destination_pos_qua_sum/cluster_size);

                    }
                    else{

                        destination_std_span= 0.0;
                        destination_std_pos=0.0;

                    }

                    score=calculate_score(cluster_size, mean(source_std_span,destination_std_span),mean(source_std_pos,destination_std_pos),mean(source_average_end-source_average_start,destination_average_end-destination_average_start),average_confidence_level, options.support_read_num);

                    if(score>options.min_SV_cluster_score){

                        consolidated_cluters.add(new SignatureClusterBilocal(SignaturePotential_Trans_Inters.getSourceContig(),(int)(Math.round(source_average_start)),(int)(Math.round(source_average_end)),SignaturePotential_Trans_Inters.getDestinationContig(),(int)(Math.round(destination_average_start)),(int)(Math.round(destination_average_end)),score,cluster.size(),cluster,SignaturePotential_Trans_Inters.getType(),mean(source_std_span,destination_std_span),mean(source_std_pos,destination_std_pos)));

                    }

                }

            }

        }

        return consolidated_cluters;

    }

    private static double mean(double v, double v1) {

        return (v+v1)/2;

    }

    private static class contig_position implements Comparator<SignatureClusterUniLocal>{

        public int compare(SignatureClusterUniLocal o1,SignatureClusterUniLocal o2){

            if(o1.contig.length()>o2.contig.length()){

                return 1;

            }
            else if(o1.contig.length()==o2.contig.length()&&o1.contig.compareTo(o2.contig)>0){

                return 1;

            }
            else if(o1.contig.length()==o2.contig.length()&&o1.contig.compareTo(o2.contig)==0&&(o1.end_position_on_ref+o1.start_position_on_ref)/2>(o2.end_position_on_ref+o2.start_position_on_ref)/2){

                return 1;

            }
            else if(o1.contig.length()==o2.contig.length()&&o1.contig.compareTo(o2.contig)==0&&(o1.end_position_on_ref+o1.start_position_on_ref)/2 == (o2.end_position_on_ref+o2.start_position_on_ref)/2){

                return 0;

            }
            else{

                return -1;

            }

        }

    }

    private static class bilocal_contig_position implements Comparator<SignatureClusterBilocal>{

        public int compare(SignatureClusterBilocal o1,SignatureClusterBilocal o2){

            if(o1.source_contig.length()>o2.source_contig.length()){

                return 1;

            }
            else if(o1.source_contig.length()==o2.source_contig.length()&&o1.source_contig.compareTo(o2.source_contig)>0){

                return 1;

            }
            else if(o1.source_contig.length()==o2.source_contig.length()&&o1.source_contig.compareTo(o2.source_contig)==0&&(o1.source_end+o1.source_start)/2>(o2.source_end+o2.source_start)/2){

                return 1;

            }
            else if(o1.source_contig.length()==o2.source_contig.length()&&o1.source_contig.compareTo(o2.source_contig)==0&&(o1.source_end+o1.source_start)/2 == (o2.source_end+o2.source_start)/2){

                return 0;

            }
            else{

                return -1;

            }

        }

    }

}
