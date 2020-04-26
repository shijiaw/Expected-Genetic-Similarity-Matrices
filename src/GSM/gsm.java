package GSM;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import fig.basic.IOUtils;
import goblin.Taxon;
import nuts.io.CSV;
import nuts.io.IO;
import pty.RootedTree;
import pty.UnrootedTree;

public class gsm {
	public static void main(String[] args) {
				
		String inputTree = "/Users/oudomame/Expected-Genetic-Similarity-Matrices/Example/tree1.newick";
		String OutputFileName = "/Users/oudomame/Expected-Genetic-Similarity-Matrices/Example/K.csv";
		PrintWriter KOut = IOUtils.openOutEasy(new File(OutputFileName));
		//PrintWriter KOut = IOUtils.openOutEasy(new File("/Users/oudomame/EGSM/Example/", OutputFileName));
		
		File dataFile = new File(inputTree);
		UnrootedTree goldut = UnrootedTree.fromRooted(RootedTree.Util.fromNewickString(IO.f2s(dataFile)));
		List<Taxon> leaves = goldut.leaves();
		String st = dataFile.getName();
		double[][] K = phylogeny.tree2K(goldut, leaves);
		
		for(int i = 0; i < leaves.size(); i++) {
			List<Double> Kind = new ArrayList<>();
			for(int ind = 0; ind < leaves.size(); ind++) {
				Kind.add(K[i][ind]);
			}
			KOut.println(CSV.body(Kind));
		}		
		
		KOut.close();
	}
}
