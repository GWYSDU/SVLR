import java.util.*;

public class T_I_C_SVClusters_Distinguish {

    public static List parse_Potential_Trans_Inters_clusters(List list, Parameter_Setting Parameter_Setting){

        List candidates=new LinkedList();
        List<SV_CandidateDeletion> candidateDeletions=new LinkedList <>();
        List<SV_CandidateInsertion> candidateInsertions=new LinkedList <>();
        List<SV_CandidateInversion> candidateInversions=new LinkedList <>();
        List<SV_CandidateDuplicationTandem> candidateDuplicationTandems =new LinkedList <>();
        List<SV_CandidateInterspersedDuplication> candidateInterspersedDuplications=new LinkedList <>();
        List<SV_CandidateTranslocation> candidateTranslocations=new LinkedList <>();

        List<SignatureClusterUniLocal> Deletions= (List <SignatureClusterUniLocal>) list.get(0);
        List<SignatureClusterUniLocal> Insertions= (List <SignatureClusterUniLocal>) list.get(1);
        List<SignatureClusterUniLocal> Inversions= (List <SignatureClusterUniLocal>) list.get(2);
        List<SignatureClusterBilocal> DuplicationTandems= (List <SignatureClusterBilocal>) list.get(3);
        List<SignatureClusterBilocal> potential_Interspersed_Trans= (List <SignatureClusterBilocal>) list.get(4);
        List<SignatureTranslocationBreakpoint> translocation_signatures_with_reversed_translocations=new LinkedList<>();
        List<SignatureTranslocationBreakpoint> reversed_translocations=new LinkedList <>();
        int translocation_signatures_size=Global.TranslocationBreakpoint_signatures.size();

        for(int m=0;m<translocation_signatures_size;m++){

            SignatureTranslocationBreakpoint signatureTranslocation=Global.TranslocationBreakpoint_signatures.get(m);

            reversed_translocations.add(new SignatureTranslocationBreakpoint(signatureTranslocation.contig2,signatureTranslocation.pos2,(signatureTranslocation.direction2.equals("rev")?"fwd":"rev"),signatureTranslocation.contig1,signatureTranslocation.pos1,(signatureTranslocation.direction1.equals("rev")?"fwd":"rev"),signatureTranslocation.signature,signatureTranslocation.confidence_level));

        }

        translocation_signatures_with_reversed_translocations.addAll(Global.TranslocationBreakpoint_signatures);
        translocation_signatures_with_reversed_translocations.addAll(reversed_translocations);

        List<SignatureTranslocationBreakpoint> translocations_fwdfwd=new LinkedList <>();
        List<SignatureTranslocationBreakpoint> translocations_revrev=new LinkedList <>();
        int translocation_signatures_with_reversed_translocations_size=translocation_signatures_with_reversed_translocations.size();

        for(int n=0;n<translocation_signatures_with_reversed_translocations_size;n++){

            SignatureTranslocationBreakpoint tra=translocation_signatures_with_reversed_translocations.get(n);

            if(tra.direction1.equals("fwd")&&tra.direction2.equals("fwd")){

                translocations_fwdfwd.add(tra);

            }
            else if(tra.direction1.equals("rev")&&tra.direction2.equals("rev")){

                translocations_revrev.add(tra);

            }

        }

        List<List<SignatureTranslocationBreakpoint>> translocation_partitions_fwdfwd=form_partitions(translocations_fwdfwd,Parameter_Setting.trans_partition_max_distance);
        List<List<SignatureTranslocationBreakpoint>> translocation_partitions_revrev=form_partitions(translocations_revrev,Parameter_Setting.trans_partition_max_distance);

        HashMap<String,List<List<SignatureTranslocationBreakpoint>>> translocation_partitions_fwdfwd_dict=new HashMap <>();
        HashMap<String,List<List<SignatureTranslocationBreakpoint>>> translocation_partitions_revrev_dict=new HashMap <>();

        for(int n=0;n<translocation_partitions_fwdfwd.size();n++){

            List<SignatureTranslocationBreakpoint> partition =translocation_partitions_fwdfwd.get(n);
            if(!translocation_partitions_fwdfwd_dict.containsKey(partition.get(0).contig1)){

                translocation_partitions_fwdfwd_dict.put(partition.get(0).contig1,new ArrayList<>());

            }

            translocation_partitions_fwdfwd_dict.get(partition.get(0).contig1).add(partition);
            translocation_partitions_fwdfwd_dict.put(partition.get(0).contig1,translocation_partitions_fwdfwd_dict.get(partition.get(0).contig1));

        }

        for(int n=0;n<translocation_partitions_revrev.size();n++){

            List<SignatureTranslocationBreakpoint> partition =translocation_partitions_revrev.get(n);
            if(!translocation_partitions_revrev_dict.containsKey(partition.get(0).contig1)){

                translocation_partitions_revrev_dict.put(partition.get(0).contig1,new ArrayList <>());

            }

            translocation_partitions_revrev_dict.get(partition.get(0).contig1).add(partition);
            translocation_partitions_revrev_dict.put(partition.get(0).contig1,translocation_partitions_revrev_dict.get(partition.get(0).contig1));

        }

        HashMap<String,ArrayList<Integer>> translocation_partitions_means_fwdfwd_dict=new HashMap <>();
        HashMap<String,ArrayList<Integer>> translocation_partitions_stds_fwdfwd_dict=new HashMap <>();

        for (Map.Entry<String,List<List<SignatureTranslocationBreakpoint>>>entry:translocation_partitions_fwdfwd_dict.entrySet()) {

            List<List<SignatureTranslocationBreakpoint>> partitions=entry.getValue();

            for(int m=0;m<partitions.size();m++){

                List<SignatureTranslocationBreakpoint> partition= partitions.get(m);
                long sum =0;

                for(int n=0;n<partition.size();n++){

                    sum+=partition.get(n).pos1;

                }

                if(!translocation_partitions_means_fwdfwd_dict.containsKey(entry.getKey())){

                    translocation_partitions_means_fwdfwd_dict.put(entry.getKey(),new ArrayList <>());

                }

                translocation_partitions_means_fwdfwd_dict.get(entry.getKey()).add((int)(Math.round((double)sum/partition.size())));
                translocation_partitions_means_fwdfwd_dict.put(entry.getKey(),translocation_partitions_means_fwdfwd_dict.get(entry.getKey()));

            }

            for(int m=0;m<partitions.size();m++){

                List<SignatureTranslocationBreakpoint> partition= partitions.get(m);
                long sum =0;

                for(int n=0;n<partition.size();n++){

                    sum+=Math.pow(Math.abs(partition.get(n).pos1-translocation_partitions_means_fwdfwd_dict.get(entry.getKey()).get(m)),2);

                }

                if(!translocation_partitions_stds_fwdfwd_dict.containsKey(entry.getKey())){

                    translocation_partitions_stds_fwdfwd_dict.put(entry.getKey(),new ArrayList <>());

                }

                translocation_partitions_stds_fwdfwd_dict.get(entry.getKey()).add((int)(Math.round(Math.sqrt((double)sum/partition.size()))));
                translocation_partitions_stds_fwdfwd_dict.put(entry.getKey(),translocation_partitions_stds_fwdfwd_dict.get(entry.getKey()));

            }

        }

        HashMap<String,ArrayList<Integer>> translocation_partitions_means_revrev_dict=new HashMap <>();
        HashMap<String,ArrayList<Integer>> translocation_partitions_stds_revrev_dict=new HashMap <>();

        for (Map.Entry<String,List<List<SignatureTranslocationBreakpoint>>>entry:translocation_partitions_revrev_dict.entrySet()) {

            List<List<SignatureTranslocationBreakpoint>> partitions=entry.getValue();

            for(int m=0;m<partitions.size();m++){

                List<SignatureTranslocationBreakpoint> partition= partitions.get(m);
                long sum =0;

                for(int n=0;n<partition.size();n++){

                    sum+=partition.get(n).pos1;

                }

                if(!translocation_partitions_means_revrev_dict.containsKey(entry.getKey())){

                    translocation_partitions_means_revrev_dict.put(entry.getKey(),new ArrayList <>());

                }

                translocation_partitions_means_revrev_dict.get(entry.getKey()).add((int)(Math.round((double)sum/partition.size())));
                translocation_partitions_means_revrev_dict.put(entry.getKey(),translocation_partitions_means_revrev_dict.get(entry.getKey()));
            }

            for(int m=0;m<partitions.size();m++){

                List<SignatureTranslocationBreakpoint> partition= partitions.get(m);
                long sum =0;

                for(int n=0;n<partition.size();n++){

                    sum+=Math.pow(Math.abs(partition.get(n).pos1-translocation_partitions_means_revrev_dict.get(entry.getKey()).get(m)),2);

                }

                if(!translocation_partitions_stds_revrev_dict.containsKey(entry.getKey())){

                    translocation_partitions_stds_revrev_dict.put(entry.getKey(),new ArrayList <>());

                }

                translocation_partitions_stds_revrev_dict.get(entry.getKey()).add((int)(Math.round(Math.sqrt((double)sum/partition.size()))));
                translocation_partitions_stds_revrev_dict.put(entry.getKey(),translocation_partitions_stds_revrev_dict.get(entry.getKey()));
            }

        }
        List returnlists=merge_translocations_at_insertions(Insertions,translocation_partitions_fwdfwd_dict,translocation_partitions_means_fwdfwd_dict,translocation_partitions_stds_fwdfwd_dict,translocation_partitions_revrev_dict,translocation_partitions_means_revrev_dict,translocation_partitions_stds_revrev_dict,Parameter_Setting);
        potential_Interspersed_Trans.addAll((Collection <? extends SignatureClusterBilocal>) returnlists.get(0)) ;

        for(int i=0;i<potential_Interspersed_Trans.size();i++){

            SignatureClusterBilocal potential_interspersed_trans=potential_Interspersed_Trans.get(i);
            String potential_interspersed_trans_source_contig=potential_interspersed_trans.source_contig;
            long potential_interspersed_trans_source_start=potential_interspersed_trans.source_start;
            long potential_interspersed_trans_source_end=potential_interspersed_trans.source_end;
            String potential_interspersed_trans_destination_contig=potential_interspersed_trans.dest_contig;
            long potential_interspersed_trans_destination_start=potential_interspersed_trans.dest_start;
            long potential_interspersed_trans_destination_end=potential_interspersed_trans.dest_end;
            boolean mark=false;

            for(int j=0;j<potential_Interspersed_Trans.size();j++){

                if(j!=i){

                    SignatureClusterBilocal potential_interspersed_trans_1=potential_Interspersed_Trans.get(j);
                    String potential_interspersed_trans_source_contig_1=potential_interspersed_trans_1.source_contig;
                    long potential_interspersed_trans_source_start_1=potential_interspersed_trans_1.source_start;
                    long potential_interspersed_trans_source_end_1=potential_interspersed_trans_1.source_end;
                    String potential_interspersed_trans_destination_contig_1=potential_interspersed_trans_1.dest_contig;
                    long potential_interspersed_trans_destination_start_1=potential_interspersed_trans_1.dest_start;
                    long potential_interspersed_trans_destination_end_1=potential_interspersed_trans_1.dest_end;

                    for(int n=0;n<potential_Interspersed_Trans.size();n++){

                        if(n!=j&&n!=i){

                            SignatureClusterBilocal potential_interspersed_trans_2=potential_Interspersed_Trans.get(n);
                            String potential_interspersed_trans_source_contig_2=potential_interspersed_trans_2.source_contig;
                            long potential_interspersed_trans_source_start_2=potential_interspersed_trans_2.source_start;
                            long potential_interspersed_trans_source_end_2=potential_interspersed_trans_2.source_end;
                            String potential_interspersed_trans_destination_contig_2=potential_interspersed_trans_2.dest_contig;
                            long potential_interspersed_trans_destination_start_2=potential_interspersed_trans_2.dest_start;
                            long potential_interspersed_trans_destination_end_2=potential_interspersed_trans_2.dest_end;

                            if(potential_interspersed_trans_source_contig.equals(potential_interspersed_trans_source_contig_1)&&potential_interspersed_trans_source_contig.equals(potential_interspersed_trans_source_contig_2)&&potential_interspersed_trans_destination_contig.equals(potential_interspersed_trans_destination_contig_1)&&potential_interspersed_trans_destination_contig.equals(potential_interspersed_trans_destination_contig_2)){

                                if(Math.abs(potential_interspersed_trans_source_start-potential_interspersed_trans_destination_end_1)<=50&&Math.abs(potential_interspersed_trans_source_end-potential_interspersed_trans_destination_start_2)<=50&&Math.abs(potential_interspersed_trans_destination_start-potential_interspersed_trans_source_end_1)<=50&&Math.abs(potential_interspersed_trans_destination_end-potential_interspersed_trans_source_start_2)<=50){

                                    mark=true;
                                    break;

                                }

                            }

                        }

                    }

                    if(mark){

                        potential_Interspersed_Trans.remove(i);
                        i--;
                        break;

                    }

                }

            }

        }

        Insertions= (List <SignatureClusterUniLocal>) returnlists.get(1);
        ArrayList<ArrayList> sources_infor=new ArrayList <>();
        int Potential_Trans_Inters_size=potential_Interspersed_Trans.size();
        int deletion_clusters_size=Deletions.size();
        for(int i=0;i<Potential_Trans_Inters_size;i++){

            ArrayList source_infor=new ArrayList();
            SignatureClusterBilocal Potential_Trans_Inters=potential_Interspersed_Trans.get(i);
            source_infor.add(Potential_Trans_Inters.source_contig);
            source_infor.add(Potential_Trans_Inters.source_start);
            source_infor.add(Potential_Trans_Inters.source_end);
            sources_infor.add(source_infor);

        }

        ArrayList<ArrayList> source_in_destination;
        boolean is_exist_source;

        for(int i=0;i<Potential_Trans_Inters_size;i++){

            source_in_destination=new ArrayList <>();
            SignatureClusterBilocal Potential_Trans_Inters=potential_Interspersed_Trans.get(i);
            int sources_infor_size=sources_infor.size();
            is_exist_source=false;

            for(int j=0;j<sources_infor_size;j++){

                ArrayList source_infor=sources_infor.get(j);

                if(Potential_Trans_Inters.dest_contig.equals(source_infor.get(0))&&Potential_Trans_Inters.dest_start<=(long)source_infor.get(1)&&(long)source_infor.get(2)<=Potential_Trans_Inters.dest_end){

                    source_in_destination.add(source_infor);
                    is_exist_source=true;

                }

            }

            if(is_exist_source){

                Collections.sort(source_in_destination,new source_in_destination_Cororder());
                long coverage_source_in_destination=get_coverage_source_in_destination(source_in_destination);
                long dis_detination=Potential_Trans_Inters.dest_end-Potential_Trans_Inters.dest_start;
                long dis_source=Potential_Trans_Inters.source_end-Potential_Trans_Inters.source_start;

                if(dis_source>=0&&dis_detination>=0&&dis_detination<100000){

                    if((dis_detination-coverage_source_in_destination-dis_source)>=20){

                        candidateTranslocations.add(new SV_CandidateTranslocation(Potential_Trans_Inters.source_contig,Potential_Trans_Inters.source_start,Potential_Trans_Inters.source_end,Potential_Trans_Inters.dest_contig,Potential_Trans_Inters.dest_start,Potential_Trans_Inters.dest_end,Potential_Trans_Inters.score));

                    }
                    else {

                        candidateInterspersedDuplications.add(new SV_CandidateInterspersedDuplication(Potential_Trans_Inters.source_contig,Potential_Trans_Inters.source_start,Potential_Trans_Inters.source_end,Potential_Trans_Inters.dest_contig,Potential_Trans_Inters.dest_start,Potential_Trans_Inters.dest_end,Potential_Trans_Inters.score,false));

                    }

                }

            }
            else{

                long dis_detination=Potential_Trans_Inters.dest_end-Potential_Trans_Inters.dest_start;
                long dis_source=Potential_Trans_Inters.source_end-Potential_Trans_Inters.source_start;

                if(dis_source>=0&&dis_detination>=0&&dis_detination<100000){

                    if((dis_detination-dis_source)>=20){

                        candidateTranslocations.add(new SV_CandidateTranslocation(Potential_Trans_Inters.source_contig,Potential_Trans_Inters.source_start,Potential_Trans_Inters.source_end,Potential_Trans_Inters.dest_contig,Potential_Trans_Inters.dest_start,Potential_Trans_Inters.dest_end,Potential_Trans_Inters.score));

                    }
                    else{

                        candidateInterspersedDuplications.add(new SV_CandidateInterspersedDuplication(Potential_Trans_Inters.source_contig,Potential_Trans_Inters.source_start,Potential_Trans_Inters.source_end,Potential_Trans_Inters.dest_contig,Potential_Trans_Inters.dest_start,Potential_Trans_Inters.dest_end,Potential_Trans_Inters.score,false));

                    }

                }

            }

        }

        for(int i=0;i<deletion_clusters_size;i++){

            SignatureClusterUniLocal Deletion=Deletions.get(i);
            candidateDeletions.add(new SV_CandidateDeletion(Deletion.contig,Deletion.start_position_on_ref,Deletion.end_position_on_ref,Deletion.score,Deletion.ave_span));

        }

        int insertion_clusters_size=Insertions.size();

        for(int i=0;i<insertion_clusters_size;i++){

            SignatureClusterUniLocal Insertion=Insertions.get(i);
            candidateInsertions.add(new SV_CandidateInsertion(Insertion.contig,Insertion.start_position_on_ref,Insertion.end_position_on_ref,Insertion.score,Insertion.ave_span));

        }

        int inversion_clusters_size=Inversions.size();

        for(int i=0;i<inversion_clusters_size;i++){

            SignatureClusterUniLocal Inversion=Inversions.get(i);
            candidateInversions.add(new SV_CandidateInversion(Inversion.contig,Inversion.start_position_on_ref,Inversion.end_position_on_ref,Inversion.score));

        }

        int duplication_tandem_clusters_size=DuplicationTandems.size();

        for(int n=0;n<duplication_tandem_clusters_size;n++){

            SignatureClusterBilocal tandem_duplication=DuplicationTandems.get(n);
            int copies=(int)(Math.round((double)(tandem_duplication.getDestinationEnd()-tandem_duplication.getDestinationStart())/(tandem_duplication.getSourceEnd()-tandem_duplication.getSourceStart())));
            candidateDuplicationTandems.add(new SV_CandidateDuplicationTandem(tandem_duplication.source_contig,tandem_duplication.source_start,tandem_duplication.source_end,copies,tandem_duplication.score));

        }

        Collections.sort(candidateInterspersedDuplications,new int_duplication_sort());
        Collections.sort(candidateDuplicationTandems,new tan_dup_sort());

        candidates.add(candidateDeletions);
        candidates.add(candidateInsertions);
        candidates.add(candidateInversions);
        candidates.add(candidateDuplicationTandems);
        candidates.add(candidateInterspersedDuplications);
        candidates.add(candidateTranslocations);

        return candidates;

    }
    private static class int_duplication_sort implements Comparator<SV_CandidateInterspersedDuplication> {

        public int compare(SV_CandidateInterspersedDuplication o1,SV_CandidateInterspersedDuplication o2){

            if(o1.contig2.length()>o2.contig2.length()){

                return 1;

            }

            else if(o1.contig2.length()==o2.contig2.length()&&o1.contig2.compareTo(o2.contig2)>0){

                return 1;

            }
            else if(o1.contig2.length()==o2.contig2.length()&&o1.contig2.compareTo(o2.contig2)==0&&o1.destination_start_position_on_ref>o2.destination_start_position_on_ref){

                return 1;

            }
            else if(o1.contig2.length()==o2.contig2.length()&&o1.contig2.compareTo(o2.contig2)==0&&o1.destination_start_position_on_ref==o2.destination_start_position_on_ref&&o1.destination_end_position_on_ref>o2.destination_end_position_on_ref){

                return 1;

            }
            else if(o1.contig2.length()==o2.contig2.length()&&o1.contig2.compareTo(o2.contig2)==0&&o1.destination_start_position_on_ref==o2.destination_start_position_on_ref&&o1.destination_end_position_on_ref==o2.destination_end_position_on_ref){

                return 0;

            }
            else{

                return  -1;

            }

        }

    }

    private static class tan_dup_sort implements Comparator<SV_CandidateDuplicationTandem> {

        public int compare(SV_CandidateDuplicationTandem o1,SV_CandidateDuplicationTandem o2){

            if(o1.contig.length()>o2.contig.length()){

                return 1;

            }

            else if(o1.contig.length()==o2.contig.length()&&o1.contig.compareTo(o2.contig)>0){

                return 1;

            }
            else if(o1.contig.length()==o2.contig.length()&&o1.contig.compareTo(o2.contig)==0&&o1.start_position_on_ref>o2.end_position_on_ref){

                return 1;

            }
            else if(o1.contig.length()==o2.contig.length()&&o1.contig.compareTo(o2.contig)==0&&o1.start_position_on_ref==o2.start_position_on_ref){

                return 0;

            }
            else{

                return  -1;

            }

        }

    }

    private static List<List<SignatureTranslocationBreakpoint>> form_partitions(List<SignatureTranslocationBreakpoint> translocations_fwdfwd, int trans_partition_max_distance) {

        List<List<SignatureTranslocationBreakpoint>> partitions=new LinkedList <>();
        List<SignatureTranslocationBreakpoint> current_partition=new LinkedList <>();
        int signaturesize=translocations_fwdfwd.size();
        Collections.sort(translocations_fwdfwd,new translocation_key_sort());

        for(int n=0;n<signaturesize;n++){

            SignatureTranslocationBreakpoint SignatureTranslocation=translocations_fwdfwd.get(n);

            if(current_partition.size()>0&&current_partition.get(0).mean_distance_to(SignatureTranslocation)>trans_partition_max_distance){

                partitions.add(current_partition);
                current_partition=new LinkedList <>();

            }

            current_partition.add(SignatureTranslocation);

        }

        if(current_partition.size()>0){

            partitions.add(current_partition);

        }

        return partitions;

    }

    private static class translocation_key_sort implements Comparator<SignatureTranslocationBreakpoint> {


        @Override
        public int compare(SignatureTranslocationBreakpoint o1, SignatureTranslocationBreakpoint o2) {

            if(o1.get_key_type().compareTo(o2.get_key_type())>0){

                return 1;

            }
            else if(o1.get_key_type().compareTo(o2.get_key_type())==0&&o1.get_key_contig().length()>o2.get_key_contig().length()){

                return 1;

            }
            else if(o1.get_key_type().compareTo(o2.get_key_type())==0&&o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())>0){

                return 1;

            }
            else if(o1.get_key_type().compareTo(o2.get_key_type())==0&&o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()>o2.get_key_position()){

                return 1;

            }
            else if(o1.get_key_type().compareTo(o2.get_key_type())==0&&o1.get_key_contig().length()==o2.get_key_contig().length()&&o1.get_key_contig().compareTo(o2.get_key_contig())==0&&o1
                    .get_key_position()==o2.get_key_position()){

                return 0;

            }
            else{

                return -1;

            }

        }

    }

    private static long get_coverage_source_in_destination(ArrayList<ArrayList> source_in_destination) {

        int source_in_destination_size=source_in_destination.size();
        long coverage=0;

        if (source_in_destination_size == 0) {

        } else if (source_in_destination_size == 1) {

            coverage=(long)source_in_destination.get(0).get(2)-(long)source_in_destination.get(0).get(1);

        } else {

            long x1s = (long)source_in_destination.get(0).get(1);
            long x1e = (long)source_in_destination.get(0).get(2);

            coverage= x1e - x1s + 1;

            for (int i = 1; i < source_in_destination_size; i++) {

                ArrayList sid = source_in_destination.get(i);
                long overlap;
                long seedlen;
                long alignlen;

                ComputeOverlap result1;
                result1 = ComputeOverlap.computeoverlap(x1s, x1e, (long)source_in_destination.get(i).get(1), (long)source_in_destination.get(i).get(2));
                overlap = result1.getOverlap();
                seedlen = result1.getSeedlen();
                alignlen = result1.getAlignlen();

                if (0 == overlap) {

                    coverage += alignlen;

                } else {

                    coverage += (alignlen - overlap);

                }
                if ((long)source_in_destination.get(i).get(1)< x1s) {

                    x1s = (long)source_in_destination.get(i).get(1);

                }
                if ((long)source_in_destination.get(i).get(2) > x1e) {

                    x1e = (long)source_in_destination.get(i).get(2);

                }

            }

        }

        return coverage;

    }

    private static class source_in_destination_Cororder implements Comparator<ArrayList> {

        public int compare(ArrayList o1, ArrayList o2) {

            if((long)o1.get(1)>(long)o2.get(1)){

                return 1;
            }
            else if((long)o1.get(1)==(long)o2.get(1)){


                return 0;
            }
            else{

                return -1;

            }

        }

    }

    private static List merge_translocations_at_insertions(List <SignatureClusterUniLocal> insertions, HashMap <String, List <List <SignatureTranslocationBreakpoint>>> translocation_partitions_fwdfwd_dict, HashMap <String, ArrayList <Integer>> translocation_partitions_means_fwdfwd_dict, HashMap <String, ArrayList <Integer>> translocation_partitions_stds_fwdfwd_dict, HashMap <String, List <List <SignatureTranslocationBreakpoint>>> translocation_partitions_revrev_dict, HashMap <String, ArrayList <Integer>> translocation_partitions_means_revrev_dict, HashMap <String, ArrayList <Integer>> translocation_partitions_stds_revrev_dict, Parameter_Setting parameter_setting) {

        List returnlists=new ArrayList();
        List<SignatureClusterBilocal> potential_interspersed_trans_clusters=new LinkedList <>();

        for(int n=0;n<insertions.size();n++){

            SignatureClusterUniLocal ins= insertions.get(n);

            String ins_contig=ins.getSourceContig();
            long ins_start=ins.getSourceStart();
            long ins_end=ins.getSourceEnd();
            int closest_to_start_fwdfwd_index;
            int closest_to_start_fwdfwd_mean;
            int closest_to_start_revrev_index;
            int closest_to_start_revrev_mean;

            try{

                closest_to_start_fwdfwd_index=get_closest_index(translocation_partitions_means_fwdfwd_dict.get(ins_contig),ins_start);
                closest_to_start_fwdfwd_mean=translocation_partitions_means_fwdfwd_dict.get(ins_contig).get(closest_to_start_fwdfwd_index);
                closest_to_start_revrev_index=get_closest_index(translocation_partitions_means_revrev_dict.get(ins_contig),ins_start);
                closest_to_start_revrev_mean=translocation_partitions_means_revrev_dict.get(ins_contig).get(closest_to_start_revrev_index);

            }
            catch (Exception e){

                e.printStackTrace();
                continue;

            }

            if(Math.abs(closest_to_start_fwdfwd_mean-ins_start)<=parameter_setting.trans_sv_max_distance&&Math.abs(closest_to_start_revrev_mean-ins_start)<=parameter_setting.trans_sv_max_distance){

                List<signature> translocation_partitions_fwdfwd=new LinkedList <>();
                List<SignatureTranslocationBreakpoint> translocation_partition_fwdfwd =translocation_partitions_fwdfwd_dict.get(ins_contig).get(closest_to_start_fwdfwd_index);
                int partitionsize=translocation_partition_fwdfwd.size();

                for(int m=0;m<partitionsize;m++){

                    SignatureTranslocationBreakpoint signature =translocation_partition_fwdfwd.get(m);
                    translocation_partitions_fwdfwd.add(new signature(signature.contig2,signature.pos2));

                }

                List<List<signature>> destinations_from_start_fwdfwd=cluster_position_simple(translocation_partitions_fwdfwd,parameter_setting.trans_destination_partition_max_distance);

                List<signature> translocation_partitions_revrev=new LinkedList <>();
                List<SignatureTranslocationBreakpoint> translocation_partition_revrev =translocation_partitions_revrev_dict.get(ins_contig).get(closest_to_start_revrev_index);
                partitionsize=translocation_partition_revrev.size();

                for(int m=0;m<partitionsize;m++){

                    SignatureTranslocationBreakpoint signature =translocation_partition_revrev.get(m);
                    translocation_partitions_revrev.add(new signature(signature.contig2,signature.pos2));

                }

                List<List<signature>> destinations_from_start_revrev=cluster_position_simple(translocation_partitions_revrev,parameter_setting.trans_destination_partition_max_distance);

                if(destinations_from_start_fwdfwd.size()==1&&destinations_from_start_revrev.size()==1){

                    long sum_fwdfwd=0;

                    for(int m=0;m<destinations_from_start_fwdfwd.get(0).size();m++){

                        sum_fwdfwd+=destinations_from_start_fwdfwd.get(0).get(m).pos;

                    }

                    signature destination_from_start_fwdfwd = new signature(destinations_from_start_fwdfwd.get(0).get(0).contig,Math.round(sum_fwdfwd/destinations_from_start_fwdfwd.get(0).size()));

                    for(int m=0;m<destinations_from_start_fwdfwd.get(0).size();m++){

                        sum_fwdfwd+=Math.pow(Math.abs(destinations_from_start_fwdfwd.get(0).get(m).pos-destination_from_start_fwdfwd.pos),2);

                    }

                    int destination_from_start_fwdfwd_std=Math.round(sum_fwdfwd/destinations_from_start_fwdfwd.get(0).size());

                    long sum_revrev=0;

                    for(int m=0;m<destinations_from_start_revrev.get(0).size();m++){

                        sum_revrev+=destinations_from_start_revrev.get(0).get(m).pos;

                    }

                    signature destination_from_start_revrev = new signature(destinations_from_start_revrev.get(0).get(0).contig,Math.round(sum_revrev/destinations_from_start_revrev.get(0).size()));

                    for(int m=0;m<destinations_from_start_revrev.get(0).size();m++){

                        sum_revrev+=Math.pow(Math.abs(destinations_from_start_revrev.get(0).get(m).pos-destination_from_start_revrev.pos),2);

                    }

                    int destination_from_start_revrev_std=Math.round(sum_revrev/destinations_from_start_revrev.get(0).size());

                    double distance= Math.abs(destination_from_start_revrev.pos-destination_from_start_fwdfwd.pos);

                    if(destination_from_start_fwdfwd.contig.equals(destination_from_start_revrev.contig)&&0.95<=(ins_end-ins_start+1)/(distance+1)&&(ins_end-ins_start+1)/(distance+1)<=1.1){

                        List members=new LinkedList();
                        members.addAll(ins.members);
                        members.addAll(translocation_partitions_fwdfwd_dict.get(ins_contig).get(closest_to_start_fwdfwd_index));
                        members.addAll(translocation_partitions_revrev_dict.get(ins_contig).get(closest_to_start_revrev_index));
                        List<Long> translocation_distances=new ArrayList <>();
                        translocation_distances.add(Math.abs(closest_to_start_fwdfwd_mean-ins_start));
                        translocation_distances.add(Math.abs(closest_to_start_revrev_mean-ins_start));
                        List<Integer> translocation_stds=new ArrayList <>();
                        translocation_stds.add(translocation_partitions_stds_fwdfwd_dict.get(ins_contig).get(closest_to_start_fwdfwd_index));
                        translocation_stds.add(translocation_partitions_stds_revrev_dict.get(ins_contig).get(closest_to_start_revrev_index));
                        List<Integer> destination_stds=new ArrayList <>();
                        destination_stds.add(destination_from_start_fwdfwd_std);
                        destination_stds.add(destination_from_start_revrev_std);
                        double score=calculate_score_insertion(ins.score,translocation_distances,translocation_stds,destination_stds);
                        potential_interspersed_trans_clusters.add(new SignatureClusterBilocal(destination_from_start_revrev.contig,Math.min(destination_from_start_revrev.pos,destination_from_start_fwdfwd.pos),Math.max(destination_from_start_revrev.pos,destination_from_start_fwdfwd.pos),ins_contig,ins_start,(long)(ins_start+distance),score,members.size(),members,"ins_dup",ins.std_span,ins.std_pos));
                        insertions.remove(n);
                        n--;

                    }

                }

            }

        }

        returnlists.add(potential_interspersed_trans_clusters);
        returnlists.add(insertions);
        return returnlists;

    }

    private static double calculate_score_insertion(double main_score, List<Long> translocation_distances, List<Integer> translocation_stds, List<Integer> destination_stds) {

        long td0=Math.max(0,100-translocation_distances.get(0))/100;
        long td1=Math.max(0,100-translocation_distances.get(1))/100;

        long ts0=Math.max(0,100-translocation_stds.get(0))/100;
        long ts1=Math.max(0,100-translocation_stds.get(1))/100;

        long ds0=Math.max(0,100-destination_stds.get(0))/100;
        long ds1=Math.max(0,100-destination_stds.get(1))/100;

        double product=(main_score/100)*td0*td1*ts0*ts1*ds0*ds1;
        double final_score=Math.pow(product,1/7)*100;

        return final_score;

    }

    private static Integer get_closest_index(ArrayList <Integer> input_list, long input_number) {

        if(input_list.size()<1){

            return null;

        }

        int pos=bisect_left(input_list,input_number);

        if(pos==0){

            return 0;

        }
        else if(pos==input_list.size()){

            return input_list.size()-1;

        }
        else {

            int before =input_list.get(pos-1);
            int after = input_list.get(pos);

            if(after-input_number<input_number-before){

                return pos;

            }
            else{

                return pos-1;

            }

        }


    }

    private static int bisect_left(ArrayList <Integer> input_list, long input_number) {

        int left=0;
        int right=input_list.size()-1;

        int middle=0;

        while(left<=right){

            middle=(right-left)/2+left;

            if(input_number==(long)input_list.get(middle)){

                return middle-1;

            }
            else if(input_number<(long)input_list.get(middle)){

                right=middle-1;

            }
            else{

                left=middle+1;

            }

        }

        if(input_number<(long)input_list.get(middle)){

            return middle;

        }
        else{

            return middle+1;

        }

    }

    private static ArrayList<List<signature>> cluster_position_simple(List <signature> translocation_partition_fwdfwd, int max_delta) {

        Collections.sort(translocation_partition_fwdfwd,new sort_position());
        ArrayList<List<signature>> partitions=new ArrayList <>();
        ArrayList<signature> current_partition=new ArrayList <>();

        for(int m=0;m<translocation_partition_fwdfwd.size();m++){

            signature position=translocation_partition_fwdfwd.get(m);

            if(current_partition.size()<1){

                current_partition.add(position);
                continue;

            }

            if(distance_positions(current_partition.get(0),position)>max_delta){

                partitions.add(current_partition);

                while(current_partition.size()>0&&distance_positions(current_partition.get(0),position)>max_delta){

                    current_partition.remove(0);

                }

            }

            current_partition.add(position);

        }

        if(current_partition.size()>0){

            partitions.add(current_partition);

        }

        return partitions;

    }

    private static class sort_position implements Comparator <signature> {
        @Override
        public int compare(signature o1, signature o2) {

            if(o1.contig.length()>o2.contig.length()){

                return 1;

            }
            else  if(o1.contig.length()==o2.contig.length()&&o1.contig.compareTo(o2.contig)>0){

                return 1;

            }

            else if(o1.contig.length()==o2.contig.length()&&o1.contig.compareTo(o2.contig)==0&&o1.pos>o2.pos){

                return 1;

            }
            else if(o1.contig.length()==o2.contig.length()&&o1.contig.compareTo(o2.contig)==0&&o1.pos==o2.pos){

                return 0;

            }
            else{

                return -1;

            }

        }

    }

    private static double distance_positions(signature signature, signature position) {

        if(!signature.contig.equals(position.contig)){

            return Double.POSITIVE_INFINITY;

        }
        else{

            return Math.abs(signature.pos-position.pos);

        }

    }

}
class signature{

    String contig;
    long pos;

    public signature(String contig,long pos){

        this.contig=contig;
        this.pos=pos;

    }

}
