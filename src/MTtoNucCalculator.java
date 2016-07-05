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

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.text.DecimalFormat;
import java.util.Iterator;

/**
 * Created by alex on 01.11.14.
 */
public class MTtoNucCalculator {
    private final SAMFileReader inputSam;
    private FileWriter fw;
    private BufferedWriter bfw;
    private long referenceLength = 0;
    private int mtlength = 0;
    private int nuclearreads = 0;
    private int mitochondrialreads = 0;
    private int lengthofmtreads = 0;
    private int lengthofnucreads = 0;
    private String mtidentifier = "";


    public MTtoNucCalculator(File f, String outputpath, String mtidentifier) throws IOException {
        inputSam = new SAMFileReader(f);
        inputSam.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
        referenceLength = inputSam.getFileHeader().getSequenceDictionary().getReferenceLength();
        this.mtidentifier = mtidentifier;
        mtlength = getMTLength(inputSam.getFileHeader().getTextHeader());

        fw = new FileWriter(new File(outputpath));
        bfw = new BufferedWriter(fw);
    }


    public static void main(String[] args) throws IOException {
        if (args.length != 3) {
            System.err.println("Please provide the (coordinate) sorted input SAM File and the output path, as well as the MT identifier. No further parameters" +
                    "are necessary! \n");
            System.err.println("Make sure that your input file has an appropriate SAM/BAM file header, or the SAMRecords will be set to '*'! \n");
            System.exit(1);
        } else {
            File f = new File(args[0]);
            String outputpath = f.getAbsolutePath() + ".mtnucratio";

            MTtoNucCalculator mTtoNucCalculator = new MTtoNucCalculator(f, outputpath, args[2]);
            mTtoNucCalculator.readSAMFile();
        }
    }

    /**
     * This Method reads a SAM File and parses the input
     * Currently, only merged Reads with the "M" Flag in front are checked for Duplicates.
     * R/F Flags are simply written into output File, also other "non-flagged" ones.
     *
     * @throws java.io.IOException
     */
    private void readSAMFile() throws IOException {
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
                    lengthofnucreads += curr.getReadLength();
                }
            }
        }
        //Write our values now into output file
        bfw.write("# of reads on mitochondrium: " + mitochondrialreads + "\n");
        Double mtcoverage = mitochondrialreads * (Double.valueOf(lengthofmtreads) / Double.valueOf(mitochondrialreads))/mtlength;
        bfw.write("AVG Coverage on MT: " + mtcoverage + "\n");
        bfw.write("# of reads on nuclear chromosomes: " + nuclearreads + "\n");
        Double nuccoverage = nuclearreads*(Double.valueOf(lengthofnucreads / Double.valueOf(nuclearreads)))/referenceLength;
        bfw.write("AVG Coverage on nuclear chromosomes: " + nuccoverage + "\n");
        DecimalFormat df = new DecimalFormat("##.##");

        if(mitochondrialreads > 0 && nuclearreads > 0 ){
            bfw.write("mt/nuc Ratio: " + df.format(Double.valueOf(mtcoverage) / Double.valueOf(nuccoverage)) + "\n");
        } else if (mitochondrialreads > 0) {
            bfw.write("mt/nuc Ratio: NF\n");
        }  else {
            bfw.write("mt/nuc Ratio: " + 0.0 + "\n");
        }
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
}
