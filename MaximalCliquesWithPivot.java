

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

class  Vertex implements Comparable<Vertex>{

    int x;
    int degree;

    public Vertex(){

        degree=0;

    }
    ArrayList<Vertex> nbrs=new ArrayList <>();

    public int getX(){

        return x;

    }

    public void setX(int x){

        this.x=x;

    }

    public int getDegree(){

        return degree;

    }

    public void setDegree(int degree){

        this.degree=degree;

    }

    public ArrayList<Vertex> getNbrs(){

        return nbrs;

    }

    public void setNbrs(ArrayList<Vertex> nbrs){

        this.nbrs=nbrs;

    }

    public void addNbr(Vertex y){

        if(!this.getNbrs().contains(y)){

            this.nbrs.add(y);

            if(!y.getNbrs().contains(this)){

                y.getNbrs().add(this);
                y.degree++;

            }

            this.degree++;

        }

    }

    public void removeNbr(Vertex y){

        this.nbrs.remove(y);

        if(y.getNbrs().contains(y)){

            y.getNbrs().remove(this);
            y.degree--;

        }

        this.degree--;

    }

    @Override
    public int compareTo(Vertex o) {

        if(this.degree<o.degree){

            return -1;

        }

        if(this.degree>o.degree){

            return 1;

        }

        return 0;

    }

}

public class MaximalCliquesWithPivot {

    int nodesCount;
    ArrayList<Vertex> graph;

    public MaximalCliquesWithPivot(int nodesCount ){

        this.nodesCount=nodesCount;
        graph=new ArrayList <>();
        for(int i=0;i<nodesCount;i++){

            Vertex V=new Vertex();
            V.setX(i);
            graph.add(V);

        }

    }


    public void Deletion_creatGraph(List <SignatureDeletion> partition_sample, Parameter_Setting options){

        try{

            if(partition_sample!=null&&partition_sample.size()==nodesCount){

                for(int i=0; i<nodesCount; i++){

                    for(int j=0; j<nodesCount; j++){

                        if(i!=j){

                            if(partition_sample.get(i).span_loc_distance(partition_sample.get(j),options.distance_normalizer)<=options.cluter_max_distance){

                                graph.get(i).addNbr(graph.get(j));
                            }

                        }

                    }

                }

            }

        }
        catch (Exception e){

            e.printStackTrace();
            throw e;

        }

    }

    public void Insertion_creatGraph(List <SignatureInsertion> partition_sample, Parameter_Setting options, long largest_indel_size){

        try{

            if(partition_sample!=null&&partition_sample.size()==nodesCount){

                for(int i=0; i<nodesCount; i++){

                    for(int j=0; j<nodesCount; j++){

                        if(i!=j){

                            if(partition_sample.get(i).span_loc_distance(partition_sample.get(j),options.distance_normalizer)<=options.cluter_max_distance){

                                graph.get(i).addNbr(graph.get(j));

                            }

                        }

                    }

                }

            }

        }
        catch (Exception e){

            e.printStackTrace();
            throw e;

        }

    }

    public void Inversion_creatGraph(List <SignatureInversion> partition_sample, Parameter_Setting options, long largest_indel_size) {

        try{

            if(partition_sample!=null&&partition_sample.size()==nodesCount){

                for(int i=0; i<nodesCount; i++){

                    for(int j=0; j<nodesCount; j++){

                        if(i!=j){

                            if(partition_sample.get(i).span_loc_distance(partition_sample.get(j),options.distance_normalizer)<=options.cluter_max_distance){

                                graph.get(i).addNbr(graph.get(j));

                            }

                        }

                    }

                }

            }

        }
        catch (Exception e){

            e.printStackTrace();
            throw e;

        }

    }

    public void Duplication_Tandem_creatGraph( List <SignatureDuplicationTandem> partition_sample, Parameter_Setting options){

        try{

            if(partition_sample!=null&&partition_sample.size()==nodesCount){

                for(int i=0; i<nodesCount; i++){

                    for(int j=0; j<nodesCount; j++){

                        if(i!=j){

                            if(partition_sample.get(i).span_loc_distance(partition_sample.get(j),options.distance_normalizer)<=options.cluter_max_distance){

                                graph.get(i).addNbr(graph.get(j));

                            }

                        }

                    }

                }

            }

        }
        catch (Exception e) {

            e.printStackTrace();
            throw e;

        }

    }

    public void potential_interspersed_trans_creatGraph(List<SignaturePotential_Trans_Inters> partition_sample, Parameter_Setting options){

        try{

            if(partition_sample!=null&&partition_sample.size()==nodesCount){

                for(int i=0; i<nodesCount; i++){

                    for(int j=0; j<nodesCount; j++){

                        if(i!=j){

                            if(partition_sample.get(i).span_loc_distance(partition_sample.get(j),options.distance_normalizer)<=options.cluter_max_distance){

                                graph.get(i).addNbr(graph.get(j));

                            }

                        }

                    }

                }

            }

        }
        catch (Exception e){

            e.printStackTrace();
            throw e;

        }

    }

    ArrayList<Vertex> getNbrs(Vertex v){

        int i=v.getX();
        return graph.get(i).nbrs;

    }

    ArrayList<Vertex> intersect(ArrayList<Vertex> arlFirst,ArrayList<Vertex> arlSecond){

        ArrayList<Vertex> arlHold =new ArrayList<Vertex> (arlFirst);
        arlHold.retainAll(arlSecond);
        return arlHold;

    }

    ArrayList<Vertex> union(ArrayList<Vertex> arlFirst,ArrayList<Vertex> arlSecond){

        ArrayList<Vertex> arlHold=new ArrayList <>(arlFirst);
        arlHold.addAll(arlSecond);
        return arlHold;

    }

    ArrayList<Vertex> removeNbrs(ArrayList<Vertex> arlFirst,Vertex V){

        ArrayList<Vertex> arlHold=new ArrayList <>(arlFirst);
        arlHold.removeAll(V.getNbrs());
        return arlHold;

    }

    void Bron_kerboshWithPivot(ArrayList<Vertex> R,ArrayList<Vertex> P,ArrayList<Vertex> X,String pre){

        if(P.size()==0&&X.size()==0){

            ArrayList<Vertex> R_copy= (ArrayList <Vertex>) R.clone();
            Global.cliques.add(R_copy);
            return ;

        }

//        while (P.size()>0){
//
//            Vertex  V=getMinDegreeVertex(P);
//            R.add(V);
//            Bron_kerboshWithPivot(R,intersect(P,getNbrs(V)),intersect(X,getNbrs(V)),pre+"\t");
//            R.remove(V);
//            P.remove(V);
//
//            ArrayList<Vertex> v_nbrs=V.getNbrs();
//            int nbrs_num=v_nbrs.size();
//
//            for(int i=0;i<nbrs_num;i++){
//
//                v_nbrs.get(i).setDegree(v_nbrs.get(i).getDegree()-1);
//
//            }
//
//            X.add(V);
//
//        }
        ArrayList<Vertex>P1=new ArrayList <>(P);
        Vertex u=getMaxDegreeVertex(union(P,X));
        P=removeNbrs(P,u);

        for(Vertex V:P){

            R.add(V);
            Bron_kerboshWithPivot(R,intersect(P1,getNbrs(V)),intersect(X,getNbrs(V)),pre+"\t");
            R.remove(V);
            P1.remove(V);
            X.add(V);
        }

    }

    void Bron_KerboschpivotExecute(){

        ArrayList<Vertex> X=new ArrayList<Vertex>();
        ArrayList<Vertex> R=new ArrayList <>();
        ArrayList<Vertex> P=new ArrayList<Vertex>(graph);
        Bron_kerboshWithPivot(R,P,X,"");

    }

    private Vertex getMaxDegreeVertex(ArrayList<Vertex> union){

        Collections.sort(union);
        return union.get(union.size()-1);

    }

    private Vertex getMinDegreeVertex(ArrayList<Vertex> union) {

//        int mindegree=union.get(0).degree;
//        Vertex V=union.get(0);
//
//        for(int i=1;i<union.size();i++){
//
//            int degree=union.get(i).degree;
//
//            if(mindegree>degree){
//
//                mindegree=degree;
//                V=union.get(i);
//
//            }
//
//        }
//
//        return V;

        Collections.sort(union);
        return union.get(0);

    }

}
