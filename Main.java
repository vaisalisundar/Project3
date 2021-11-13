import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Main {

    public static void main(String[] args) throws NullPointerException {

        List<Float> pData = new ArrayList<>();
        List<Integer> noOfBumps = new ArrayList<>();
        List<Integer> cellsProcessed = new ArrayList<>();
        List<Integer> finalLength = new ArrayList<>();
        List<Boolean> solvabilty = new ArrayList<>();
        List<Integer> trajectoryLen = new ArrayList<>();
        List<Float> trajectoryDivideLength = new ArrayList<>();
        for (int l = 0; l < 1; l++) {
            int dim = 101;
            int values[][] = new int[dim][dim];
            int dummyValues[][] = new int[dim][dim];
            double beliefState[][] = new double[dim][dim];


            double p = 0.3;


            //creating a dummy grid with no blocked cells
            for (int rn = 0; rn < dummyValues.length; rn++) {
                for (int cl = 0; cl < dummyValues[rn].length; cl++) {
                    dummyValues[rn][cl] = 0;
                    System.out.print(dummyValues[rn][cl]);
                }
                System.out.println();
            }
            System.out.println("Done printing the dummy grid");


            //terrain grid
            for (int i = 0; i < values.length; i++) {
                // do the for in the row according to the column size
                for (int j = 0; j < values[i].length; j++) {
                    float random = (float) Math.random();
                    //System.out.println(random);

                    if (random >= p) {
                        if(random>=0.3 && random<=0.54){
                            //flat terrain
                            values[i][j] = 2;
                        }else if(random>=0.54 && random<=0.77){
                            //hilly terrain
                            values[i][j] = 3;
                        }else if(random>=0.77 && random<=1.0){
                            //forest terain
                            values[i][j] = 4;
                        }

                    } else {
                        values[i][j] = 1;
                    }


                    System.out.print(values[i][j]);

                }
                // add a new line
                System.out.println();
            }

            //Assign agent and target to cells
            int min=1;
            int max = dim-1;
//            int x1=5;
//            int y1=4;
//            int x2=3;
//            int y2=4;
            int x1=-1;
            int y1=-1;
            int x2=-1;
            int y2=-1;
            for(int i=0;i<2;i++){
                if(i==0){
                    do{
                        x1 = min + (int)(Math.random() * ((max - min) + 1));
                        y1 = min + (int)(Math.random() * ((max - min) + 1));
                    }while(values[x1][y1]==1);
                }


                if(i==1){
                    do{
                        x2 = min + (int)(Math.random() * ((max - min) + 1));
                        y2 = min + (int)(Math.random() * ((max - min) + 1));
                    }while(values[x2][y2]==1);
                }

            }
            System.out.println("source : "+x1+","+y1);
            System.out.println("target : "+x2+","+y2);
            Solution sol = new Solution();
            //grid is the actual maze
             int[][] grid = values;
//            int[][] grid ={{4,2,4,4,3,1,1,4,4,3},
//                    {4,4,4,2,1,4,3,1,4,3},
//                    {4,3,3,4,1,2,1,1,4,3},
//                    {1,4,2,1,4,2,2,2,3,2},
//                    {1,2,1,2,3,2,1,4,3,4},
//                    {2,4,2,4,2,1,3,1,1,1},
//                    {4,3,1,2,1,3,3,3,1,1},
//                    {2,2,4,1,3,4,4,4,1,1},
//                    {3,4,3,3,2,3,2,1,1,3},
//                    {3,3,2,4,1,2,1,2,1,2}};
            //int[][] grid ={{0,0,1,1,1},{1,0,1,1,1},{0,0,1,1,1},{0,1,1,1,1},{0,0,0,0,0}};


            //newGrid is the maze with no blockage
            int[][] newGrid = dummyValues;


            int length = sol.shortestPathBinaryMatrix(grid,x1,y1,x2,y2);
            System.out.println(length);
            System.out.println();

            if (length != -1) {
                // float ans = sol.RFAstarTest(newGrid,0,0,1,grid);
                 sol.AgentSix(newGrid, x1, y1, 1, grid,beliefState,dim,x2,y2);
                System.out.println("source : "+x1+","+y1);
                System.out.println("target : "+x2+","+y2);
               // System.out.println(ans);
                //pData.add(p);
                //float result = (float)ans/(float)length;
                // System.out.println(ans);
                // trajectoryLen.add(ans);
                int cells = sol.noOfCellsProcessed;
                //noOfBumps.add(ans);
               // finalLength.add(ans);
                //trajectoryLen.add(ans);
                cellsProcessed.add(cells);
                // trajectoryDivideLength.add(ans);
            }


            //System.out.println(cells);

            //No.of.cells processed in RFA
            //int cellsRFA = sol.cellsProcessedRFA()-1;

            //pData.add(p);
            //lengthData.add(length);

            //cellsProcessedRFA.add(cellsRFA);


        }
        FileWriter locFile1 = null;
        FileWriter locFile2 = null;
        //     FileWriter locFile3 = null;
//        FileWriter locFile4 = null;
        //       FileWriter locFile5 = null;
        // FileWriter locFile6 = null;


        try{
            locFile1 = new FileWriter("finalLen.txt");
//            locFile2 = new FileWriter("lengthData.txt");
            locFile2 = new FileWriter("cellsProcessed.txt");
//            locFile4 = new FileWriter("solvability.txt");
//            locFile5 = new FileWriter("cellProcessedRFA.txt");
            // locFile6 = new FileWriter("trajectoryDivideLen.txt");


            for(int i = 0;i<pData.size();i++){
                locFile1.write(finalLength.get(i)+"\n");
                locFile2.write(cellsProcessed.get(i)+"\n");
                //      locFile3.write(cellsProcessed.get(i)+"\n");
//                locFile4.write(solvabilty.get(i)+"\n");
//                locFile5.write(cellsProcessedRFA.get(i)+"\n");
                //   locFile6.write(trajectoryDivideLength.get(i)+"\n");
            }
            //  locFile6.close();
//            locFile5.close();
//            locFile4.close();
//            locFile3.close();
            locFile2.close();
            locFile1.close();
        }catch (IOException e){
            System.out.println("In catch block");
            e.printStackTrace();
        }


    }



}


class Solution {
    double flatFNR = 0.2;
    double hillyFNR = 0.5;
    double forestFNR = 0.8;
    Map<List<Integer>, List<Integer>> pathMapRFA = new HashMap<List<Integer>, List<Integer>>();
    List<Candidate> track = new ArrayList<Candidate>();


    int noOfCellsProcessed = 0;

    boolean[][] visitedRFA = new boolean[8][8];

    //Queue<Candidate> pqRFA = new PriorityQueue<>((a, b) -> a.totalEstimate - b.totalEstimate);
    // Candidate represents a possible option for travelling to the cell
    // at (row, col).
    class Candidate {

        public int row;
        public int col;
        public int distanceSoFar;
        public int totalEstimate;


        public Candidate(int row, int col, int distanceSoFar, int totalEstimate) {
            this.row = row;
            this.col = col;
            this.distanceSoFar = distanceSoFar;
            this.totalEstimate = totalEstimate;
        }


    }

    class Equation {

        public List<Candidate> candidateList;
        public int c;



        public Equation(List<Candidate> candidateList , int c) {
            this.candidateList = candidateList;
            this.c = c;
            //inferenceKnowledge.put(c,candidateList);


        }


    }


    private static final int[][] directions =
            new int[][]{{-1, 0}, {0, -1}, {0, 1}, {1, 0}};

    //Covering diagonal neighbors -Agent-3
    private static final int[][] directionsAgent3 =
            new int[][]{{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};



    public int shortestPathBinaryMatrix(int[][] grid,int x1,int y1,int x2,int y2) {

        // Firstly, we need to check that the start and target cells are open.
        if (grid[x1][y1] == 1 || grid[x2][y2] == 1) {
            return -1;
        }

        // Set up the A* search.
        Queue<Candidate> pq = new PriorityQueue<>((a, b) -> a.totalEstimate - b.totalEstimate);
        pq.add(new Candidate(x1, y1, 1, estimateAgent6(x1, y1, grid,x2,y2)));

        boolean[][] visited = new boolean[grid.length][grid[0].length];
        Map<List<Integer>, List<Integer>> pathMap = new HashMap<List<Integer>, List<Integer>>();


        // Carry out the A* search.
        while (!pq.isEmpty()) {

            Candidate best = pq.remove();
            // int counter = cellsProcessed();


            // Is this for the target cell?
            if (best.row == x2 && best.col == y2) {

                List<Integer> goalnode = new ArrayList<>();
                goalnode.add(best.row);
                goalnode.add(best.col);
                System.out.print(goalnode);
                System.out.print("->");
                List<Integer> result = pathMap.get(goalnode);
                for (int i = 0; i < best.distanceSoFar; i++) {
                    System.out.print(result);
                    System.out.print("->");
                    result = pathMap.get(result);

                }

                return best.distanceSoFar;
            }

            // Have we already found the best path to this cell?
            if (visited[best.row][best.col]) {
                continue;
            }

            visited[best.row][best.col] = true;
            // pathQueue.add({new int[]{best.row,best.col},{});

            for (int[] neighbour : getNeighbours(best.row, best.col, grid)) {
                int neighbourRow = neighbour[0];
                int neighbourCol = neighbour[1];

                // This check isn't necessary for correctness, but it greatly
                // improves performance.
                if (visited[neighbourRow][neighbourCol]) {
                    continue;
                }

                // Otherwise, we need to create a Candidate object for the option
                // of going to neighbor through the current cell.
                int newDistance = best.distanceSoFar + 1;
                int totalEstimate = newDistance + estimateAgent6(neighbourRow, neighbourCol, grid,x2,y2);
                Candidate candidate =
                        new Candidate(neighbourRow, neighbourCol, newDistance, totalEstimate);
                pq.add(candidate);
                List<Integer> neigh = new ArrayList<>();
                neigh.add(neighbourRow);
                neigh.add(neighbourCol);
                List<Integer> parent = new ArrayList<>();
                parent.add(best.row);
                parent.add(best.col);
                if (pathMap.containsKey(neigh)) {
                    continue;
                } else {
                    pathMap.put(neigh, parent);
                }


            }
        }
        // The target was unreachable.
        return -1;
    }

    private List<int[]> getNeighbours(int row, int col, int[][] grid) {
        List<int[]> neighbours = new ArrayList<>();
        for (int i = 0; i < directions.length; i++) {
            int newRow = row + directions[i][0];
            int newCol = col + directions[i][1];
            if (newRow < 0 || newCol < 0 || newRow >= grid.length
                    || newCol >= grid[0].length
                    || grid[newRow][newCol] == 1) {
                continue;
            }
            neighbours.add(new int[]{newRow, newCol});
        }
        return neighbours;
    }

//    // Get the best case estimate of how much further it is to the bottom-right cell.
//    private int estimate(int row, int col, int[][] grid) {
//        //Manhattan distance
//        int heuristic = Math.abs(row - (grid.length - 1)) + Math.abs(col - (grid[0].length - 1));
//        return heuristic;
//    }

    private int estimateAgent6(int row, int col, int[][] grid,int x,int y) {
        //Manhattan distance
        int heuristic = Math.abs(row - x) + Math.abs(col - y);
        return heuristic;
    }

    private void intialBielfStateSetting(double[][] beliefState,int dim){

        for(int i=0;i<beliefState.length;i++){

            for(int j=0;j<beliefState[i].length;j++){


                 beliefState[i][j] = (double)1/(dim*dim);
            }
        }
    }


    //examine function
    private int examine(int x1,int y1,int x2,int y2,double[][] beliefState,
                        int[][] trueGrid,int[][] dummyGrid,int currRow,
                        int currCol){
        //At time t + 1 you attempt to enter (x, y) and find it is blocked?
        // Px,y (t+1) = 0
        // Pi,j(t+1)=Pi,j(t)/ i=0 to dim & j=0 to dim Pi,j(t)  (Summation)
        if(trueGrid[currRow][currCol]==1){
            dummyGrid[currRow][currCol]=1;
            beliefState[currRow][currCol]=0;
            divide(beliefState,currRow,currCol);
            return 1;
        }

        if(currRow!=x2 || currCol!=y2) {
            //flat terrain
            //At time t + 1 you examine cell (x, y) of terrain type flat, and fail to find the target?
            //Px,y(t+1) = Px,y(t)* FN(Tr)
            //Pi,j(t+1) = Pi,j(t)/i=0 to dim & j=0 to dim Pi,j(t) (Summation)
            if (trueGrid[currRow][currCol] == 2) {
                examineTerrainType((double) flatFNR, beliefState, currRow, currCol);
                divide(beliefState, currRow, currCol);
                return 0;
            }
            //hilly terrain
            //At time t + 1 you examine cell (x, y) of terrain type hilly, and fail to find the target?
            //Px,y(t+1) = Px,y(t)* FN(Tr)
            //Pi,j(t+1) = Pi,j(t)/i=0 to dim & j=0 to dim Pi,j(t) (Summation)
            if (trueGrid[currRow][currCol] == 3) {
                examineTerrainType((double)hillyFNR, beliefState, currRow, currCol);
                divide(beliefState, currRow, currCol);
                return 0;
            }
            //forest Terrain
            //At time t + 1 you examine cell (x, y) of terrain type forest, and fail to find the target?
            //Px,y(t+1) = Px,y(t)* FN(Tr) /i=0 dimj=0 dimPi,j(t)
            //Pi,j(t+1) = Pi,j(t)/i=0 to dim & j=0 to dim Pi,j(t) (Summation)
            if (trueGrid[currRow][currCol] == 4) {
                examineTerrainType((double) forestFNR, beliefState, currRow, currCol);
                divide(beliefState, currRow, currCol);
                return 0;
            }
        }else{
            //target is in cell
            beliefState[currRow][currCol] = 1;
            return -1;
        }
        return 0;
    }


    //Px,y(t+1) = Px,y(t)* FN(Tr)
    private void examineTerrainType(double FNR,double[][] beliefState,int currRow,int currCol){
        double variable;
        variable = (beliefState[currRow][currCol]*FNR);
        beliefState[currRow][currCol]=variable;
    }
    //P(i,j)T/Summation
    private void divide(double[][] beliefState,int currRow,int currCol){
        double variable;
        double summation  = summation(beliefState,currRow,currCol);
        for (int rn = 0; rn < beliefState.length; rn++) {
            for (int cl = 0; cl < beliefState[rn].length; cl++) {
                variable = beliefState[rn][cl]/summation;
                beliefState[rn][cl]=variable;
            }
        }
    }

    //summation
    private double summation(double[][] beliefState,int currRow,int currCol){
        double sum =0;

        for (int rn = 0; rn < beliefState.length; rn++) {
            for (int cl = 0; cl < beliefState[rn].length; cl++) {
                sum += beliefState[rn][cl];


            }
        }
        return sum;
    }

    private double findMaximumProbability(double[][] beliefState){
        double max = Integer.MIN_VALUE;
        for(int i =0;i<beliefState.length;i++){
            for(int j=0;j<beliefState[i].length;j++){
                if(beliefState[i][j]>max){
                    max = beliefState[i][j];
                }
            }
        }
        return max;
    }

    private Candidate choosingCellWithMaxProbability(double[][] beliefState, double max, int x1, int y1, int[][] grid){
        ArrayList<Candidate> cellsWithMaxProb = new ArrayList<Candidate>();
        for(int i=0;i< beliefState.length;i++){
            for(int j=0;j< beliefState[i].length;j++){
                if(beliefState[i][j]==max){
//                    if(i==x1 && j==y1){
//                        continue;
//                    }
                    cellsWithMaxProb.add(new Candidate(i,j,estimateAgent6(i,j,grid,x1,y1),0));
                }
            }
        }

        if(cellsWithMaxProb.size()!=1){
            int min = Integer.MAX_VALUE;
            for(int i=0;i<cellsWithMaxProb.size();i++){
                if(cellsWithMaxProb.get(i).distanceSoFar<min ){
                    min = cellsWithMaxProb.get(i).distanceSoFar;
                }
            }

            for(int j=0;j<cellsWithMaxProb.size();j++){
                if(cellsWithMaxProb.get(j).distanceSoFar > min){
                    cellsWithMaxProb.remove(j);
                    j = -1;
                }
            }
            int lowLimit=0;
            int upperLimit = cellsWithMaxProb.size()-1;
            if(cellsWithMaxProb.size()!=1){
                int randomNumber = lowLimit + (int)(Math.random() * ((upperLimit - lowLimit) + 1));
                return cellsWithMaxProb.get(randomNumber);
            }
        }


        return cellsWithMaxProb.get(0);

    }

    private Candidate choosingCellWithMaxProbabilitySPL(double[][] beliefState, double max, int x1, int y1, int[][] grid,Candidate cellTobeIgnored){
        ArrayList<Candidate> cellsWithMaxProb = new ArrayList<Candidate>();
        for(int i=0;i< beliefState.length;i++){
            for(int j=0;j< beliefState[i].length;j++){
                if(beliefState[i][j]==max){
                    if(i==cellTobeIgnored.row && j==cellTobeIgnored.col){
                        continue;
                    }
                    cellsWithMaxProb.add(new Candidate(i,j,estimateAgent6(i,j,grid,x1,y1),0));
                }
            }
        }

        if(cellsWithMaxProb.size()!=1){
            int min = Integer.MAX_VALUE;
            for(int i=0;i<cellsWithMaxProb.size();i++){
                if(cellsWithMaxProb.get(i).distanceSoFar<min ){
                    min = cellsWithMaxProb.get(i).distanceSoFar;
                }
            }

            for(int j=0;j<cellsWithMaxProb.size();j++){
                if(cellsWithMaxProb.get(j).distanceSoFar > min){
                    cellsWithMaxProb.remove(j);
                    j = -1;
                }
            }
            int lowLimit=0;
            int upperLimit = cellsWithMaxProb.size()-1;
            if(cellsWithMaxProb.size()!=1){
                int randomNumber = lowLimit + (int)(Math.random() * ((upperLimit - lowLimit) + 1));
                return cellsWithMaxProb.get(randomNumber);
            }
        }


        return cellsWithMaxProb.get(0);

    }



    public static int cellsProcessed = 0;


    public static int cellsProcessedRFA = 0;



    public int RFAstar(int[][] newGrid, int row, int col, int distanceSoFar,int x2,int y2) {
        int counter = 0;
        int resultCount =0;

        // Firstly, we need to check that the start and target cells are open.
        if (newGrid[row][col] == 1 || newGrid[x2][y2] == 1) {
            return -1;
        }

        // Set up the A* search.
        Queue<Candidate> pq = new PriorityQueue<>((a, b) -> a.totalEstimate - b.totalEstimate);
        pq.add(new Candidate(row, col, distanceSoFar, estimateAgent6(0, 0, newGrid,x2,y2)));

        boolean[][] visited = new boolean[newGrid.length][newGrid[0].length];
        Map<Candidate, Candidate> pathMap = new HashMap<Candidate, Candidate>();


        // Carry out the A* search.
        while (!pq.isEmpty()) {

            Candidate best = pq.remove();
            noOfCellsProcessed++;
            counter++;
            resultCount++;

            // Is this for the target cell?
            if (best.row == x2 && best.col == y2) {


                if (track.isEmpty()) {
                    Candidate result = pathMap.get(best);
                    System.out.print(best.row+","+best.col);
                    System.out.print("->");
                    track.add(best);
                    if(resultCount>1) {
                        Candidate source = new Candidate(row, col, distanceSoFar, estimateAgent6(row, col, newGrid, x2, y2));
                        // for (int i = 0; i < best.distanceSoFar - 1; i++) {
                        while (result.row != source.row || result.col != source.col) {
                            System.out.print(result.row + "," + result.col);
                            track.add(result);
                            System.out.print("->");
                            result = pathMap.get(result);

                        }
                        System.out.println(source.row + "," + source.col);
                        track.add(source);
                    }
                } else {
                    track.clear();
                    Candidate result = pathMap.get(best);
                    System.out.print(best.row + "," + best.col);
                    System.out.print("->");
                    track.add(best);
                    Candidate source = new Candidate(row, col, distanceSoFar, estimateAgent6(row, col, newGrid, x2, y2));
                    //             for (int i = 0; i <7; i++) {
                    if (resultCount > 1) {
                        while (result.row != source.row || result.col != source.col) {
                            System.out.print(result.row + "," + result.col);
                            track.add(result);
                            System.out.print("->");
                            result = pathMap.get(result);

                        }
                        track.add(source);
                        //System.out.println(source.row+","+ source.col);
                        System.out.println();

                    }
                }
                return best.distanceSoFar;
                //return counter;
            }

            // Have we already found the best path to this cell?
            if (visited[best.row][best.col]) {
                continue;
            }

            visited[best.row][best.col] = true;
            // pathQueue.add({new int[]{best.row,best.col},{});

            for (int[] neighbour : getNeighbours(best.row, best.col, newGrid)) {
                int neighbourRow = neighbour[0];
                int neighbourCol = neighbour[1];

                // This check isn't necessary for correctness, but it greatly
                // improves performance.
                if (visited[neighbourRow][neighbourCol]) {
                    continue;
                }

                // Otherwise, we need to create a Candidate object for the option
                // of going to neighbor through the current cell.
                int newDistance = best.distanceSoFar + 1;
                int totalEstimate = newDistance + estimateAgent6(neighbourRow, neighbourCol, newGrid,x2,y2);
                Candidate candidate =
                        new Candidate(neighbourRow, neighbourCol, newDistance, totalEstimate);
                pq.add(candidate);

                if (!pathMap.containsKey(candidate)) {
                    pathMap.put(candidate, best);
                }
                // pathMap.put(candidate, best);



            }
        }
        // The target was unreachable.
        return -1;
    }

    //AgentSix - blindfolded Probabilistic sensing
    public boolean AgentSix(int[][] newGrid, int x1, int y1, int distanceSoFar, int[][] grid,double[][] beliefState,int dim,int x2,int y2) {
        int noOfBumps = 0;
        int totalCounter=0;
        int trajectoryLen = 0;
        boolean[][] visited = new boolean[newGrid.length][newGrid[0].length];
        intialBielfStateSetting(beliefState,dim);
        double max = findMaximumProbability(beliefState);
        Candidate cellWithMaxProb = choosingCellWithMaxProbability(beliefState,max,x1,y1,grid);
        int ans = RFAstar(newGrid, x1, y1, distanceSoFar,cellWithMaxProb.row,cellWithMaxProb.col);
        totalCounter =  ans+totalCounter;
        Collections.reverse(track);

        for (int i = 0; i < track.size(); i++) {
            int examineResult;
            Candidate canValue = track.get(i);
            if(visited[canValue.row][canValue.col]!=true && i==(track.size()-1)){
                 examineResult = examine(x1,y1,x2,y2,beliefState,grid,newGrid,canValue.row, canValue.col);
                visited[canValue.row][canValue.col]=true;
                if(examineResult==1){
                    double newMax = findMaximumProbability(beliefState);
                    if(max<newMax) {
                        max = newMax;
                        cellWithMaxProb = choosingCellWithMaxProbability(beliefState, newMax, x1, y1, grid);
                    }
                    Candidate newSource = track.get(i - 1);
                    int answer = RFAstar(newGrid, newSource.row, newSource.col, 0,cellWithMaxProb.row,cellWithMaxProb.col);
                    while(answer==-1){
                        beliefState[cellWithMaxProb.row][cellWithMaxProb.col]=0;
                        divide(beliefState, cellWithMaxProb.row, cellWithMaxProb.col);
                        newMax = findMaximumProbability(beliefState);
                        if(max<newMax) {
                            max = newMax;
                            //Candidate copy = cellWithMaxProb;
                            cellWithMaxProb = choosingCellWithMaxProbability(beliefState, newMax, x1, y1, grid);
                        }
                        answer = RFAstar(newGrid, newSource.row, newSource.col, 0,cellWithMaxProb.row,cellWithMaxProb.col);
                    }
                    Collections.reverse(track);
                    i = -1;
                }else if(examineResult==0){
                    double newMax = findMaximumProbability(beliefState);
                    if(max<newMax){
                        max = newMax;
                        //Candidate copy = cellWithMaxProb;
                        cellWithMaxProb = choosingCellWithMaxProbability(beliefState, newMax, x1, y1, grid);

                        int answer = RFAstar(newGrid, x1, y1, 0,cellWithMaxProb.row,cellWithMaxProb.col);
                        while(answer==-1){
                            beliefState[cellWithMaxProb.row][cellWithMaxProb.col]=0;
                            divide(beliefState, cellWithMaxProb.row, cellWithMaxProb.col);
                            newMax = findMaximumProbability(beliefState);
                            if(max<newMax) {
                                max = newMax;
                                //Candidate copy = cellWithMaxProb;
                                cellWithMaxProb = choosingCellWithMaxProbability(beliefState, newMax, x1, y1, grid);
                            }
                             answer = RFAstar(newGrid, x1, y1, 0,cellWithMaxProb.row,cellWithMaxProb.col);
                        }
                        Collections.reverse(track);
                        i = -1;


                    }
                } else if(examineResult==-1){
                    System.out.println("target found");
                    // RFAstar(newGrid, x1, y1, 1,cellWithMaxProb.row,cellWithMaxProb.col);
                    return true;
                }
            }else{
                if(grid[canValue.row][canValue.col]==1){
                    newGrid[canValue.row][canValue.col]=1;
                    beliefState[canValue.row][canValue.col]=0;
                    divide(beliefState, canValue.row, canValue.col);
                    Candidate newSource = track.get(i - 1);
                    double newMax = findMaximumProbability(beliefState);
                    if(max<newMax) {
                        max = newMax;
                        cellWithMaxProb = choosingCellWithMaxProbability(beliefState, newMax, x1, y1, grid);
                    }

                    int answer = RFAstar(newGrid, newSource.row, newSource.col, 0,cellWithMaxProb.row,cellWithMaxProb.col);
                    while(answer==-1){
                        beliefState[cellWithMaxProb.row][cellWithMaxProb.col]=0;
                        divide(beliefState, cellWithMaxProb.row, cellWithMaxProb.col);
                         newMax = findMaximumProbability(beliefState);
                        if(max<newMax) {
                            max = newMax;
                            //Candidate copy = cellWithMaxProb;
                            cellWithMaxProb = choosingCellWithMaxProbability(beliefState, newMax, x1, y1, grid);
                        }
                        answer = RFAstar(newGrid, newSource.row, newSource.col, 0,cellWithMaxProb.row,cellWithMaxProb.col);
                    }
                    Collections.reverse(track);
                    i = -1;
//                    int answer = RFAstar(newGrid, newSource.row, newSource.col, 0,cellWithMaxProb.row,cellWithMaxProb.col);
////                    if(answer==-1){
////
////                            cellWithMaxProb = choosingCellWithMaxProbabilitySPL(beliefState,max, x1, y1, grid,cellWithMaxProb);
////
////                        RFAstar(newGrid, x1, y1, 0,cellWithMaxProb.row,cellWithMaxProb.col);
////                   }
//                    Collections.reverse(track);
//                    i = -1;

                }
            }




        }
        System.out.println("target not found");
        return false;

    }


}