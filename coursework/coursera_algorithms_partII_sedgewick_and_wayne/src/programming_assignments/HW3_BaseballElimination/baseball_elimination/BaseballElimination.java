package programming_assignments.HW3_BaseballElimination.baseball_elimination;

import edu.princeton.cs.algs4.*;
import java.util.*;

public class BaseballElimination {
  // create a baseball division from given filename in format specified below
  private int n;
  private Map<String, Integer> team2no;
  private String[] no2team;
  private int[] wins;
  private int[] losses;
  private int[] remaining;
  private int[][] scheduled;
  private boolean[] isEliminated;
  private List[] certificate;

  /**
   * read input file
   * figure out which team is eliminated
   *
   * input format example
   * # of team
   * team #win #loss #remaining_game (# # # #)remaining_game_against_each_team
   4
   Atlanta       83 71  8  0 1 6 1
   Philadelphia  80 79  3  1 0 0 2
   New_York      78 78  6  6 0 0 0
   Montreal      77 82  3  1 2 0 0
   *
   * note the number of remaining games is not necessarily equal
   * to the total number of remaining games against each team
   * because a team may play with teams outside division
   *
   *
   * @param filename
   */
  public BaseballElimination(String filename) {
    In in = new In(filename);
    n = in.readInt();
    if (n < 1) throw new IllegalArgumentException("n must be >=1 ");
    team2no = new HashMap<>();
    no2team = new String[n];
    wins = new int[n];
    losses = new int[n];
    remaining = new int[n];
    scheduled = new int[n][n];
    isEliminated = new boolean[n];
    certificate = new List[n];

    for (int i = 0; i < n; i++) {
      String team = in.readString();
      team2no.put(team, i);
      no2team[i] = team;
      wins[i] = in.readInt();
      losses[i] = in.readInt();
      remaining[i] = in.readInt();
      for (int j = 0; j < n; j++) {
        scheduled[i][j] = in.readInt();
      }
    }

    /*
    for ith team, construct flow network without edges to i
    set source -> game link with capacity as remaining games
    set game -> team capacity inf
    set team -> sink as # of wins needed to tie with team i

    source(0) ---> (game vertices) ---> (team vertices) ----> sink

    vertice indices:
    source: 0
    team j: j+1
    game between j and k: sequential
    sink: N - 1

     */
    for (int i = 0; i < n; i++) {
      //every time we initialize a flow network with all teams in it
      FlowNetwork net = new FlowNetwork(1 + (n - 1) * (n - 2) / 2 + n + 1);
      //set up team -> sink
      for (int j = 0; j < n; j++) {
        if (j == i)
          continue;
        int minWinsForTie = wins[i] + remaining[i] - wins[j];
        if (minWinsForTie < 0) {
          isEliminated[i] = true;
          certificate[i] = new ArrayList<String>();
          certificate[i].add(no2team[j]);
          break;
        }
        net.addEdge(new FlowEdge(j + 1, net.V() - 1, wins[i] + remaining[i] - wins[j]));
      }
      if (isEliminated[i])
        continue;

      int nodeCount = n + 1;
      //set up source -> game -> team
      for (int j = 0; j < n; j++) {
        if (j == i)
          continue;
        for (int k = j + 1; k < n; k++) {
          if (k == i)
            continue;
          net.addEdge(new FlowEdge(0, nodeCount, scheduled[j][k]));
          net.addEdge(new FlowEdge(nodeCount, j + 1, Double.POSITIVE_INFINITY));
          net.addEdge(new FlowEdge(nodeCount, k + 1, Double.POSITIVE_INFINITY));
          nodeCount++;
        }
      }
      FordFulkerson maxFlow = new FordFulkerson(net, 0, net.V() - 1);
      for (FlowEdge e : net.adj(0)) {
        if (e.residualCapacityTo(e.other(0)) > 0) {
          //there are unplayed games
          isEliminated[i] = true;
          certificate[i] = new ArrayList<String>();
          //now figure out certificate of elimination
          for (int j = 0; j < n; j++) {
            if (j == i)
              continue;
            if (maxFlow.inCut(j + 1))
              certificate[i].add(no2team[j]);
          }
          break;
        }
      }
    }
  }

  /**
   number of teams
    */
  public int numberOfTeams() {
    return n;
  }

  private int getTeamIndex(String team) {
    if (!team2no.containsKey(team)) throw new IllegalArgumentException("no such team");
    return team2no.get(team);
  }

  /**
   all teams
    */
  public Iterable<String> teams() {
    return team2no.keySet();
  }

  /**
   number of wins for given team
    */
  public int wins(String team) {
    return wins[getTeamIndex(team)];
  }

  /**
   number of losses for given team
    */
  public int losses(String team) {
    return losses[getTeamIndex(team)];
  }

  /**
   number of remaining games for given team
    */
  public int remaining(String team) {
    return remaining[getTeamIndex(team)];
  }

  /**
   number of remaining games between team1 and team2
    */
  public int against(String team1, String team2) {
    return scheduled[getTeamIndex(team1)][getTeamIndex(team2)];
  }

  /** is given team eliminated?
   *
   * @param team
   * @return
   */
  public boolean isEliminated(String team) {
    return isEliminated[getTeamIndex(team)];
  }

  /** subset R of teams that eliminates given team; null if not eliminated
   *
   * @param team
   * @return
   */
  public Iterable<String> certificateOfElimination(String team) {
    return certificate[getTeamIndex(team)];
  }

  public static void main(String[] args) {
    BaseballElimination division = new BaseballElimination(args[0]);
    for (String team : division.teams()) {
      if (division.isEliminated(team)) {
        StdOut.print(team + " is eliminated by the subset R = { ");
        for (String t : division.certificateOfElimination(team)) {
          StdOut.print(t + " ");
        }
        StdOut.println("}");
      } else {
        StdOut.println(team + " is not eliminated");
      }
    }
  }
}