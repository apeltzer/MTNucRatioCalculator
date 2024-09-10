/*
 * Copyright (c) 2016. MTNucRatioCalculator Alexander Peltzer
 * This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import htsjdk.samtools.*;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Iterator;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

/**
 * Created by alex on 01.11.14.
 */
public class MTtoNucCalculator {
    private final SamReader inputSam;
    private FileWriter fw;
    private FileWriter jsonfw;
    private BufferedWriter bfw;
    private long referenceLength = 0;
    private int mtlength = 0;
    private int nuclearreads = 0;
    private int mitochondrialreads = 0;
    private double lengthofmtreads = 0;
    private double lengthofnucreads = 0;
    private String mtidentifier = "";
    private HashMap<String, Object> json_map = new HashMap<>();
    private static final String VERSION = "0.7.1";


    public MTtoNucCalculator(File f, String outputpath, String mtidentifier) throws IOException {
        inputSam = SamReaderFactory.make().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX).validationStringency(ValidationStringency.LENIENT).open(f);
        referenceLength = inputSam.getFileHeader().getSequenceDictionary().getReferenceLength();
        this.mtidentifier = mtidentifier;
        mtlength = getMTLength(inputSam.getFileHeader().getTextHeader());
	
        fw = new FileWriter(new File(outputpath));
        jsonfw = new FileWriter(new File(outputpath+"mtnuc.json"));
        bfw = new BufferedWriter(fw);
    }


    public static void main(String[] args) throws IOException {
        if (args.length != 2) {
            System.err.println("Version: " + VERSION);
            System.err.println("Please provide the (coordinate) sorted input SAM File, as well as the MT identifier. No further parameters" +
                    " are necessary! \n");
            System.err.println("Make sure that your input file has an appropriate SAM/BAM file header, or the SAMRecords will be set to '*'! \n");
            System.exit(1);
        } else {
            File f = new File(args[0]);
            String outputpath = f.getAbsolutePath() + ".mtnucratio";

            MTtoNucCalculator mTtoNucCalculator = new MTtoNucCalculator(f, outputpath, args[1]);
            mTtoNucCalculator.readSAMFile(f.getAbsolutePath());
        }
    }

    /**
     * This Method reads a SAM File and parses the input
     * Currently, only merged Reads with the "M" Flag in front are checked for Duplicates.
     * R/F Flags are simply written into output File, also other "non-flagged" ones.
     *
     * @throws java.io.IOException
     */
    private void readSAMFile(String filename) throws IOException {
        Iterator it = inputSam.iterator();
        while (it.hasNext()) {
            SAMRecord curr = (SAMRecord) it.next();
            //If is mapped, then count nuclear or mitochondrial read into statistics!
            if (!(curr.getReadUnmappedFlag())) {
                if (curr.getReferenceName().contains(this.mtidentifier)) {
                    mitochondrialreads++;
                    lengthofmtreads += curr.getReadLength();

                } else {
                    nuclearreads++;
                    lengthofnucreads += curr.getReadLength(); //TODO this is a problematic thing...
                }
            }
        }
        //Write our values now into output file
        HashMap<String, Object> metrics_map = new HashMap<String, Object>();
        bfw.write("# of reads on mitochondrium: " + mitochondrialreads + "\n");
        metrics_map.put("mtreads", mitochondrialreads);
        Double mtcoverage = mitochondrialreads * (Double.valueOf(lengthofmtreads) / Double.valueOf(mitochondrialreads))/mtlength;
        bfw.write("AVG Coverage on MT: " + mtcoverage + "\n");
        metrics_map.put("mt_cov_avg", mtcoverage);
        bfw.write("# of reads on nuclear chromosomes: " + nuclearreads + "\n");
        metrics_map.put("nucreads", nuclearreads);
        Double nuccoverage = nuclearreads*(Double.valueOf(lengthofnucreads / Double.valueOf(nuclearreads)))/referenceLength;
        bfw.write("AVG Coverage on nuclear chromosomes: " + nuccoverage + "\n");
        metrics_map.put("nuc_cov_avg", nuccoverage);
        DecimalFormat df = new DecimalFormat("##.##");

        if(mitochondrialreads > 0 && nuclearreads > 0 ){
            String ratio = df.format(Double.valueOf(mtcoverage) / Double.valueOf(nuccoverage));
            bfw.write("mt/nuc Ratio: " + ratio + "\n");
            metrics_map.put("mt_nuc_ratio", ratio);
        } else if (mitochondrialreads > 0) {
            bfw.write("mt/nuc Ratio: NF\n");
            metrics_map.put("mt_nuc_ratio", "NF");
        }  else {
            bfw.write("mt/nuc Ratio: " + 0.0 + "\n");
            metrics_map.put("mt_nuc_ratio", 0.0);
        }

        json_map.put("metrics", metrics_map);

        writeJSON(filename);
        bfw.flush();
        bfw.close();
    }

    /**
     * Extracts the exact SAMFileheader from the genome, to be more generic
     * @param samTextHeader
     * @return
     */
    private int getMTLength(String samTextHeader){
        String[] splitAt = samTextHeader.split("@"); //@SQ Lines are split now
        for(String s : splitAt){
            if(s.contains(this.mtidentifier)){
                //Then we found our MT ;-)
                String[] splitAgain = s.split("LN:");
                int out = Integer.valueOf(splitAgain[1].trim());
                return out;
            }
        }
        return -1;
    }


    /**
     * writes all generated output statistics to YAML output file
     * YAML has several key,value pairs (HashMaps) that all contain information created in DamageProfiler
     * sample_name = the name of the sample
     *
     *
     * @throws IOException
     */
    public void writeJSON(String input) throws IOException {
        //Add Sample Name to yaml
        String sampleName = input.split("/")[input.split("/").length-1];


        //Add Metadata to JSON output
        HashMap<String, Object> meta_map = new HashMap<>();

        meta_map.put("sample_name",sampleName);
        meta_map.put("tool_name", "mtnuccalculator");
        meta_map.put("version", VERSION);

        json_map.put("metadata", meta_map);

        //Now add the stuff we need to have too

        Gson gson = new GsonBuilder().setPrettyPrinting().serializeSpecialFloatingPointValues().create();
        String json = gson.toJson(json_map);

        jsonfw.write(json);
        jsonfw.flush();
        jsonfw.close();
    }

}
