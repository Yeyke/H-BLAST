/* $Id: blast_args.cpp 383536 2012-12-14 21:12:01Z rafanovi $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/*****************************************************************************

File name: blast_args.cpp

Author: Jason Papadopoulos

******************************************************************************/

/** @file blast_args.cpp
 * convert blast-related command line
 * arguments into blast options
*/

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = "$Id: blast_args.cpp 383536 2012-12-14 21:12:01Z rafanovi $";
#endif

#include <ncbi_pch.hpp>
#include <algo/blast/api/version.hpp>
#include <algo/blast/blastinput/blast_args.hpp>
#include <algo/blast/api/blast_exception.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/api/objmgr_query_data.hpp> /* for CObjMgrQueryFactory */
#include <algo/blast/core/blast_nalookup.h>
#include <algo/blast/core/hspfilter_besthit.h>
#include <objects/scoremat/PssmWithParameters.hpp>
#include <util/format_guess.hpp>
#include <objtools/blast/seqdb_reader/seqdb.hpp>
#include <algo/blast/blastinput/blast_input.hpp>    // for CInputException
#include <connect/ncbi_connutil.h>

#include <algo/blast/api/msa_pssm_input.hpp>    // for CPsiBlastInputClustalW
#include <algo/blast/api/pssm_engine.hpp>       // for CPssmEngine

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(blast)
USING_SCOPE(objects);
USING_SCOPE(align_format);

void
IBlastCmdLineArgs::ExtractAlgorithmOptions(const CArgs& /* cmd_line_args */,
        CBlastOptions& /* options */)
{}

CProgramDescriptionArgs::CProgramDescriptionArgs(const string& program_name,
        const string& program_desc)
    : m_ProgName(program_name), m_ProgDesc(program_desc)
{}

void
CProgramDescriptionArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    // program description
    arg_desc.SetUsageContext(m_ProgName, m_ProgDesc + " " +
                             CBlastVersion().Print());
}

CTaskCmdLineArgs::CTaskCmdLineArgs(const set<string>& supported_tasks,
                                   const string& default_task)
    : m_SupportedTasks(supported_tasks), m_DefaultTask(default_task)
{
    _ASSERT( !m_SupportedTasks.empty() );
    if ( !m_DefaultTask.empty() )
    {
        _ASSERT(m_SupportedTasks.find(m_DefaultTask) != m_SupportedTasks.end());
    }
}

void
CTaskCmdLineArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("General search options");
    if ( !m_DefaultTask.empty() )
    {
        arg_desc.AddDefaultKey(kTask, "task_name", "Task to execute",
                               CArgDescriptions::eString, m_DefaultTask);
    }
    else
    {
        arg_desc.AddKey(kTask, "task_name", "Task to execute",
                        CArgDescriptions::eString);
    }
    arg_desc.SetConstraint(kTask, new CArgAllowStringSet(m_SupportedTasks));
    arg_desc.SetCurrentGroup("");

}

void
CTaskCmdLineArgs::ExtractAlgorithmOptions(const CArgs& /* cmd_line_args */,
        CBlastOptions& /* options */)
{
    // N.B.: handling of tasks occurs at the application level to ensure that
    // only relevant tasks are added (@sa CBlastnAppArgs)
}

void
CGenericSearchArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("General search options");

    // evalue cutoff
    if (!m_IsIgBlast)
    {
        arg_desc.AddDefaultKey(kArgEvalue, "evalue",
                               "Expectation value (E) threshold for saving hits ",
                               CArgDescriptions::eDouble,
                               NStr::DoubleToString(BLAST_EXPECT_VALUE));
    }
    else if (m_QueryIsProtein)
    {
        arg_desc.AddDefaultKey(kArgEvalue, "evalue",
                               "Expectation value (E) threshold for saving hits ",
                               CArgDescriptions::eDouble,
                               NStr::DoubleToString(1.0));
    }
    else
    {
        arg_desc.AddDefaultKey(kArgEvalue, "evalue",
                               "Expectation value (E) threshold for saving hits ",
                               CArgDescriptions::eDouble,
                               NStr::DoubleToString(1e-15));
    }

    // word size
    // Default values: blastn=11, megablast=28, others=3
    const string description = m_QueryIsProtein
                               ? "Word size for wordfinder algorithm"
                               : "Word size for wordfinder algorithm (length of best perfect match)";
    arg_desc.AddOptionalKey(kArgWordSize, "int_value", description,
                            CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgWordSize, m_QueryIsProtein
                           ? new CArgAllowValuesGreaterThanOrEqual(2)
                           : new CArgAllowValuesGreaterThanOrEqual(4));

    if ( !m_IsRpsBlast && !m_IsTblastx)
    {
        // gap open penalty
        arg_desc.AddOptionalKey(kArgGapOpen, "open_penalty",
                                "Cost to open a gap",
                                CArgDescriptions::eInteger);

        // gap extend penalty
        arg_desc.AddOptionalKey(kArgGapExtend, "extend_penalty",
                                "Cost to extend a gap",
                                CArgDescriptions::eInteger);
    }


    if (m_ShowPercentIdentity)
    {
        arg_desc.SetCurrentGroup("Restrict search or results");
        arg_desc.AddOptionalKey(kArgPercentIdentity, "float_value",
                                "Percent identity",
                                CArgDescriptions::eDouble);
        arg_desc.SetConstraint(kArgPercentIdentity,
                               new CArgAllow_Doubles(0.0, 100.0));
    }

    arg_desc.SetCurrentGroup("Extension options");
    // ungapped X-drop
    // Default values: blastn=20, megablast=10, others=7
    arg_desc.AddOptionalKey(kArgUngappedXDropoff, "float_value",
                            "X-dropoff value (in bits) for ungapped extensions",
                            CArgDescriptions::eDouble);

    // Tblastx is ungapped only.
    if (!m_IsTblastx)
    {
        // initial gapped X-drop
        // Default values: blastn=30, megablast=20, tblastx=0, others=15
        arg_desc.AddOptionalKey(kArgGappedXDropoff, "float_value",
                                "X-dropoff value (in bits) for preliminary gapped extensions",
                                CArgDescriptions::eDouble);

        // final gapped X-drop
        // Default values: blastn/megablast=50, tblastx=0, others=25
        arg_desc.AddOptionalKey(kArgFinalGappedXDropoff, "float_value",
                                "X-dropoff value (in bits) for final gapped alignment",
                                CArgDescriptions::eDouble);
    }

    arg_desc.SetCurrentGroup("Statistical options");
    // effective search space
    // Default value is the real size
    arg_desc.AddOptionalKey(kArgEffSearchSpace, "int_value",
                            "Effective length of the search space",
                            CArgDescriptions::eInt8);
    arg_desc.SetConstraint(kArgEffSearchSpace,
                           new CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddDefaultKey(kArgMaxHSPsPerSubject, "int_value",
                           "Override maximum number of HSPs per subject to save for ungapped searches (0 means do not override)",
                           CArgDescriptions::eInteger,
                           NStr::IntToString(kDfltArgMaxHSPsPerSubject));
    arg_desc.SetConstraint(kArgMaxHSPsPerSubject,
                           new CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.SetCurrentGroup("");
}

void
CGenericSearchArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opt)
{
    if (args[kArgEvalue])
    {
        opt.SetEvalueThreshold(args[kArgEvalue].AsDouble());
    }

    int gap_open=0, gap_extend=0;
    if (args.Exist(kArgMatrixName) && args[kArgMatrixName])
        BLAST_GetProteinGapExistenceExtendParams
        (args[kArgMatrixName].AsString().c_str(), &gap_open, &gap_extend);

    if (args.Exist(kArgGapOpen) && args[kArgGapOpen])
    {
        opt.SetGapOpeningCost(args[kArgGapOpen].AsInteger());
    }
    else if (args.Exist(kArgMatrixName) && args[kArgMatrixName])
    {
        opt.SetGapOpeningCost(gap_open);
    }

    if (args.Exist(kArgGapExtend) && args[kArgGapExtend])
    {
        opt.SetGapExtensionCost(args[kArgGapExtend].AsInteger());
    }
    else if (args.Exist(kArgMatrixName) && args[kArgMatrixName])
    {
        opt.SetGapExtensionCost(gap_extend);
    }

    if (args[kArgUngappedXDropoff])
    {
        opt.SetXDropoff(args[kArgUngappedXDropoff].AsDouble());
    }

    if (args.Exist(kArgGappedXDropoff) && args[kArgGappedXDropoff])
    {
        opt.SetGapXDropoff(args[kArgGappedXDropoff].AsDouble());
    }

    if (args.Exist(kArgFinalGappedXDropoff) && args[kArgFinalGappedXDropoff])
    {
        opt.SetGapXDropoffFinal(args[kArgFinalGappedXDropoff].AsDouble());
    }

    if (args[kArgWordSize])
    {
        if (m_QueryIsProtein && args[kArgWordSize].AsInteger() > 5)
            opt.SetLookupTableType(eCompressedAaLookupTable);
        opt.SetWordSize(args[kArgWordSize].AsInteger());
    }

    if (args[kArgEffSearchSpace])
    {
        CNcbiEnvironment env;
        env.Set("OLD_FSC", "true");
        opt.SetEffectiveSearchSpace(args[kArgEffSearchSpace].AsInt8());
    }

    if (args.Exist(kArgPercentIdentity) && args[kArgPercentIdentity])
    {
        opt.SetPercentIdentity(args[kArgPercentIdentity].AsDouble());
    }

    if (args[kArgMaxHSPsPerSubject])
    {
        const int value = args[kArgMaxHSPsPerSubject].AsInteger();
        if (value != kDfltArgMaxHSPsPerSubject)
        {
            opt.SetMaxNumHspPerSequence(value);
        }
    }
}

void
CFilteringArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("Query filtering options");

    if (m_QueryIsProtein)
    {
        arg_desc.AddDefaultKey(kArgSegFiltering, "SEG_options",
                               "Filter query sequence with SEG "
                               "(Format: '" + kDfltArgApplyFiltering + "', " +
                               "'window locut hicut', or '" + kDfltArgNoFiltering +
                               "' to disable)",
                               CArgDescriptions::eString, m_FilterByDefault
                               ? kDfltArgSegFiltering : kDfltArgNoFiltering);
        arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
                               "Apply filtering locations as soft masks",
                               CArgDescriptions::eBoolean,
                               kDfltArgLookupTableMaskingOnlyProt);
    }
    else
    {
        arg_desc.AddDefaultKey(kArgDustFiltering, "DUST_options",
                               "Filter query sequence with DUST "
                               "(Format: '" + kDfltArgApplyFiltering + "', " +
                               "'level window linker', or '" + kDfltArgNoFiltering +
                               "' to disable)",
                               CArgDescriptions::eString, m_FilterByDefault
                               ? kDfltArgDustFiltering : kDfltArgNoFiltering);
        arg_desc.AddOptionalKey(kArgFilteringDb, "filtering_database",
                                "BLAST database containing filtering elements (i.e.: repeats)",
                                CArgDescriptions::eString);

        arg_desc.AddOptionalKey(kArgWindowMaskerTaxId, "window_masker_taxid",
                                "Enable WindowMasker filtering using a Taxonomic ID",
                                CArgDescriptions::eInteger);

        arg_desc.AddOptionalKey(kArgWindowMaskerDatabase, "window_masker_db",
                                "Enable WindowMasker filtering using this repeats database.",
                                CArgDescriptions::eString);

        arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
                               "Apply filtering locations as soft masks",
                               CArgDescriptions::eBoolean,
                               kDfltArgLookupTableMaskingOnlyNucl);
    }

    arg_desc.SetCurrentGroup("");
}

void
CFilteringArgs::x_TokenizeFilteringArgs(const string& filtering_args,
                                        vector<string>& output) const
{
    output.clear();
    NStr::Tokenize(filtering_args, " ", output);
    if (output.size() != 3)
    {
        NCBI_THROW(CInputException, eInvalidInput,
                   "Invalid number of arguments to filtering option");
    }
}

void
CFilteringArgs::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& opt)
{
    if (args[kArgLookupTableMaskingOnly])
    {
        opt.SetMaskAtHash(args[kArgLookupTableMaskingOnly].AsBoolean());
    }

    vector<string> tokens;

    try
    {
        if (m_QueryIsProtein && args[kArgSegFiltering])
        {
            const string& seg_opts = args[kArgSegFiltering].AsString();
            if (seg_opts == kDfltArgNoFiltering)
            {
                opt.SetSegFiltering(false);
            }
            else if (seg_opts == kDfltArgApplyFiltering)
            {
                opt.SetSegFiltering(true);
            }
            else
            {
                x_TokenizeFilteringArgs(seg_opts, tokens);
                opt.SetSegFilteringWindow(NStr::StringToInt(tokens[0]));
                opt.SetSegFilteringLocut(NStr::StringToDouble(tokens[1]));
                opt.SetSegFilteringHicut(NStr::StringToDouble(tokens[2]));
            }
        }

        if ( !m_QueryIsProtein && args[kArgDustFiltering])
        {
            const string& dust_opts = args[kArgDustFiltering].AsString();
            if (dust_opts == kDfltArgNoFiltering)
            {
                opt.SetDustFiltering(false);
            }
            else if (dust_opts == kDfltArgApplyFiltering)
            {
                opt.SetDustFiltering(true);
            }
            else
            {
                x_TokenizeFilteringArgs(dust_opts, tokens);
                opt.SetDustFilteringLevel(NStr::StringToInt(tokens[0]));
                opt.SetDustFilteringWindow(NStr::StringToInt(tokens[1]));
                opt.SetDustFilteringLinker(NStr::StringToInt(tokens[2]));
            }
        }
    }
    catch (const CStringException& e)
    {
        if (e.GetErrCode() == CStringException::eConvert)
        {
            NCBI_THROW(CInputException, eInvalidInput,
                       "Invalid input for filtering parameters");
        }
    }

    int filter_dbs = 0;

    if (args.Exist(kArgFilteringDb) && args[kArgFilteringDb])
    {
        opt.SetRepeatFilteringDB(args[kArgFilteringDb].AsString().c_str());
        filter_dbs++;
    }

    if (args.Exist(kArgWindowMaskerTaxId) &&
            args[kArgWindowMaskerTaxId])
    {

        opt.SetWindowMaskerTaxId
        (args[kArgWindowMaskerTaxId].AsInteger());

        filter_dbs++;
    }

    if (args.Exist(kArgWindowMaskerDatabase) &&
            args[kArgWindowMaskerDatabase])
    {

        opt.SetWindowMaskerDatabase
        (args[kArgWindowMaskerDatabase].AsString().c_str());

        filter_dbs++;
    }

    if (filter_dbs > 1)
    {
        string msg =
            string("Please specify at most one of ") + kArgFilteringDb + ", " +
            kArgWindowMaskerTaxId + ", or " + kArgWindowMaskerDatabase + ".";

        NCBI_THROW(CInputException, eInvalidInput, msg);
    }
}

void
CWindowSizeArg::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("Extension options");
    // 2-hit wordfinder window size
    arg_desc.AddOptionalKey(kArgWindowSize, "int_value",
                            "Multiple hits window size, use 0 to specify "
                            "1-hit algorithm",
                            CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgWindowSize,
                           new CArgAllowValuesGreaterThanOrEqual(0));
    arg_desc.SetCurrentGroup("");
}

void
CWindowSizeArg::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& opt)
{
    if (args[kArgWindowSize])
    {
        opt.SetWindowSize(args[kArgWindowSize].AsInteger());
    }
    else
    {
        int window = -1;
        BLAST_GetSuggestedWindowSize(opt.GetProgramType(),
                                     opt.GetMatrixName(),
                                     &window);
        if (window != -1)
        {
            opt.SetWindowSize(window);
        }
    }
}

void
COffDiagonalRangeArg::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("Extension options");
    // 2-hit wordfinder off diagonal range
    arg_desc.AddDefaultKey(kArgOffDiagonalRange, "int_value",
                           "Number of off-diagonals to search for the 2nd hit, "
                           "use 0 to turn off",
                           CArgDescriptions::eInteger,
                           NStr::IntToString(kDfltOffDiagonalRange));
    arg_desc.SetConstraint(kArgOffDiagonalRange,
                           new CArgAllowValuesGreaterThanOrEqual(0));
    arg_desc.SetCurrentGroup("");
}

void
COffDiagonalRangeArg::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& opt)
{
    if (args[kArgOffDiagonalRange])
    {
        opt.SetOffDiagonalRange(args[kArgOffDiagonalRange].AsInteger());
    }
    else
    {
        opt.SetOffDiagonalRange(0);
    }
}

// Options specific to rmblastn -RMH-
void
CRMBlastNArg::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("General search options");

    arg_desc.AddDefaultKey(kArgMatrixName, "matrix_name",
                           "Scoring matrix name",
                           CArgDescriptions::eString,
                           string(""));

    arg_desc.AddFlag(kArgComplexityAdj,
                     "Use complexity adjusted scoring",
                     true);


    arg_desc.AddDefaultKey(kArgMaskLevel, "int_value",
                           "Masklevel - percentage overlap allowed per "
                           "query domain [0-101]",
                           CArgDescriptions::eInteger,
                           kDfltArgMaskLevel);
    arg_desc.SetConstraint(kArgMaskLevel,
                           new CArgAllowValuesLessThanOrEqual(101));

    arg_desc.SetCurrentGroup("");
}

// Options specific to rmblastn -RMH-
void
CRMBlastNArg::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& opt)
{
    if (args[kArgMatrixName])
    {
        opt.SetMatrixName(args[kArgMatrixName].AsString().c_str());
    }

    opt.SetComplexityAdjMode( args[kArgComplexityAdj] );

    if (args[kArgMaskLevel])
    {
        opt.SetMaskLevel(args[kArgMaskLevel].AsInteger());
    }

    if (args[kArgMinRawGappedScore])
    {
        opt.SetCutoffScore(args[kArgMinRawGappedScore].AsInteger());
    }
    else if (args[kArgUngappedXDropoff])
    {
        opt.SetCutoffScore(args[kArgUngappedXDropoff].AsInteger());
    }
}

void
CWordThresholdArg::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("General search options");
    // lookup table word score threshold
    arg_desc.AddOptionalKey(kArgWordScoreThreshold, "float_value",
                            "Minimum word score such that the word is added to the "
                            "BLAST lookup table",
                            CArgDescriptions::eDouble);
    arg_desc.SetConstraint(kArgWordScoreThreshold,
                           new CArgAllowValuesGreaterThanOrEqual(0));
    arg_desc.SetCurrentGroup("");
}

void
CWordThresholdArg::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opt)
{
    if (args[kArgWordScoreThreshold])
    {
        opt.SetWordThreshold(args[kArgWordScoreThreshold].AsDouble());
    }
    else
    {
        double threshold = -1;
        BLAST_GetSuggestedThreshold(opt.GetProgramType(),
                                    opt.GetMatrixName(),
                                    &threshold);
        if (threshold != -1)
        {
            opt.SetWordThreshold(threshold);
        }
    }
}

void
CMatrixNameArg::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("General search options");
    arg_desc.AddOptionalKey(kArgMatrixName, "matrix_name",
                            "Scoring matrix name (normally BLOSUM62)",
                            CArgDescriptions::eString);
    arg_desc.SetCurrentGroup("");
}

void
CMatrixNameArg::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& opt)
{
    if (args[kArgMatrixName])
    {
        opt.SetMatrixName(args[kArgMatrixName].AsString().c_str());
    }
}

void
CNuclArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    // TLM arg_desc.SetCurrentGroup("Nucleotide scoring options");

    arg_desc.SetCurrentGroup("General search options");
    // blastn mismatch penalty
    arg_desc.AddOptionalKey(kArgMismatch, "penalty",
                            "Penalty for a nucleotide mismatch",
                            CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMismatch,
                           new CArgAllowValuesLessThanOrEqual(0));

    // blastn match reward
    arg_desc.AddOptionalKey(kArgMatch, "reward",
                            "Reward for a nucleotide match",
                            CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMatch,
                           new CArgAllowValuesGreaterThanOrEqual(0));


    arg_desc.SetCurrentGroup("Extension options");
    arg_desc.AddFlag(kArgNoGreedyExtension,
                     "Use non-greedy dynamic programming extension",
                     true);

    arg_desc.SetCurrentGroup("");
}

void
CNuclArgs::ExtractAlgorithmOptions(const CArgs& cmd_line_args,
                                   CBlastOptions& options)
{
    if (cmd_line_args[kArgMismatch])
    {
        options.SetMismatchPenalty(cmd_line_args[kArgMismatch].AsInteger());
    }
    if (cmd_line_args[kArgMatch])
    {
        options.SetMatchReward(cmd_line_args[kArgMatch].AsInteger());
    }

    if (cmd_line_args[kArgNoGreedyExtension])
    {
        options.SetGapExtnAlgorithm(eDynProgScoreOnly);
        options.SetGapTracebackAlgorithm(eDynProgTbck);
    }
}

const string CDiscontiguousMegablastArgs::kTemplType_Coding("coding");
const string CDiscontiguousMegablastArgs::kTemplType_Optimal("optimal");
const string
CDiscontiguousMegablastArgs::kTemplType_CodingAndOptimal("coding_and_optimal");

void
CDiscontiguousMegablastArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("Extension options");
    // FIXME: this can be applied to any program, but since it was only offered
    // in megablast, we're putting it here
    arg_desc.AddOptionalKey(kArgMinRawGappedScore, "int_value",
                            "Minimum raw gapped score to keep an alignment "
                            "in the preliminary gapped and traceback stages",
                            CArgDescriptions::eInteger);

    arg_desc.SetCurrentGroup("Discontiguous MegaBLAST options");

    arg_desc.AddOptionalKey(kArgDMBTemplateType, "type",
                            "Discontiguous MegaBLAST template type",
                            CArgDescriptions::eString);
    arg_desc.SetConstraint(kArgDMBTemplateType, &(*new CArgAllow_Strings,
                           kTemplType_Coding,
                           kTemplType_Optimal,
                           kTemplType_CodingAndOptimal));
    arg_desc.SetDependency(kArgDMBTemplateType,
                           CArgDescriptions::eRequires,
                           kArgDMBTemplateLength);

    arg_desc.AddOptionalKey(kArgDMBTemplateLength, "int_value",
                            "Discontiguous MegaBLAST template length",
                            CArgDescriptions::eInteger);
    set<int> allowed_values;
    allowed_values.insert(16);
    allowed_values.insert(18);
    allowed_values.insert(21);
    arg_desc.SetConstraint(kArgDMBTemplateLength,
                           new CArgAllowIntegerSet(allowed_values));
    arg_desc.SetDependency(kArgDMBTemplateLength,
                           CArgDescriptions::eRequires,
                           kArgDMBTemplateType);

    arg_desc.SetCurrentGroup("");
}

void
CDiscontiguousMegablastArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& options)
{
    if (args[kArgMinRawGappedScore])
    {
        options.SetCutoffScore(args[kArgMinRawGappedScore].AsInteger());
    }

    if (args[kArgDMBTemplateType])
    {
        const string& type = args[kArgDMBTemplateType].AsString();
        EDiscWordType temp_type = eMBWordCoding;

        if (type == kTemplType_Coding)
        {
            temp_type = eMBWordCoding;
        }
        else if (type == kTemplType_Optimal)
        {
            temp_type = eMBWordOptimal;
        }
        else if (type == kTemplType_CodingAndOptimal)
        {
            temp_type = eMBWordTwoTemplates;
        }
        else
        {
            abort();
        }
        options.SetMBTemplateType(static_cast<unsigned char>(temp_type));
    }

    if (args[kArgDMBTemplateLength])
    {
        unsigned char tlen =
            static_cast<unsigned char>(args[kArgDMBTemplateLength].AsInteger());
        options.SetMBTemplateLength(tlen);
    }

    // FIXME: should the window size be adjusted if this is set?
}

void
CCompositionBasedStatsArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("General search options");
    // composition based statistics, keep in sync with ECompoAdjustModes
    // documentation in composition_constants.h

    string zero_opt = !m_ZeroOptDescr.empty() ?
                      (string)"    0 or F or f: " + m_ZeroOptDescr + "\n" :
                      "    0 or F or f: No composition-based statistics\n";

    string one_opt_insrt = m_Is2and3Supported ? "" : " or T or t";

    string more_opts = m_Is2and3Supported ?
                       "    2 or T or t : Composition-based score adjustment as in "
                       "Bioinformatics 21:902-911,\n"
                       "    2005, conditioned on sequence properties\n"
                       "    3: Composition-based score adjustment as in "
                       "Bioinformatics 21:902-911,\n"
                       "    2005, unconditionally\n" : "";

    string legend = (string)"Use composition-based statistics:\n"
                    "    D or d: default (equivalent to " + m_DefaultOpt + " )\n"
                    + zero_opt
                    + "    1" + one_opt_insrt + ": Composition-based statistics "
                    "as in NAR 29:2994-3005, 2001\n"
                    + more_opts;

    arg_desc.AddDefaultKey(kArgCompBasedStats, "compo", legend,
                           CArgDescriptions::eString, m_DefaultOpt);


    arg_desc.SetCurrentGroup("Miscellaneous options");
    // Use Smith-Waterman algorithm in traceback stage
    // FIXME: available only for gapped blastp/tblastn, and with
    // composition-based statistics
    arg_desc.AddFlag(kArgUseSWTraceback,
                     "Compute locally optimal Smith-Waterman alignments?",
                     true);
    arg_desc.SetCurrentGroup("");
}

/**
 * @brief Auxiliary function to set the composition based statistics and smith
 * waterman options
 *
 * @param opt BLAST options object [in|out]
 * @param comp_stat_string command line value for composition based statistics
 * [in]
 * @param smith_waterman_value command line value for determining the use of
 * the smith-waterman algorithm [in]
 * @param ungapped pointer to the value which determines whether the search
 * should be ungapped or not. It is NULL if ungapped searches are not
 * applicable
 * @param is_deltablast is program deltablast [in]
 */
static void
s_SetCompositionBasedStats(CBlastOptions& opt,
                           const string& comp_stat_string,
                           bool smith_waterman_value,
                           bool* ungapped)
{
    const EProgram program = opt.GetProgram();
    if (program == eBlastp || program == eTblastn ||
            program == ePSIBlast || program == ePSITblastn ||
            program == eRPSBlast ||
            program == eBlastx  ||  program == eDeltaBlast)
    {

        ECompoAdjustModes compo_mode = eNoCompositionBasedStats;

        switch (comp_stat_string[0])
        {
        case '0':
        case 'F':
        case 'f':
            compo_mode = eNoCompositionBasedStats;
            break;
        case '1':
            compo_mode = eCompositionBasedStats;
            break;
        case 'D':
        case 'd':
            if (program == eRPSBlast)
            {
                compo_mode = eNoCompositionBasedStats;
            }
            else if (program == eDeltaBlast)
            {
                compo_mode = eCompositionBasedStats;
            }
            else
            {
                compo_mode = eCompositionMatrixAdjust;
            }
            break;
        case '2':
            compo_mode = eCompositionMatrixAdjust;
            break;
        case '3':
            compo_mode = eCompoForceFullMatrixAdjust;
            break;
        case 'T':
        case 't':
            compo_mode = (program == eRPSBlast || program == eDeltaBlast) ?
                         eCompositionBasedStats : eCompositionMatrixAdjust;
            break;
        }

        if(program == ePSITblastn)
        {
            compo_mode = eNoCompositionBasedStats;
        }

        if (ungapped && *ungapped && compo_mode != eNoCompositionBasedStats)
        {
            NCBI_THROW(CInputException, eInvalidInput,
                       "Composition-adjusted searched are not supported with "
                       "an ungapped search, please add -comp_based_stats F or "
                       "do a gapped search");
        }

        opt.SetCompositionBasedStats(compo_mode);
        if (program == eBlastp &&
                compo_mode != eNoCompositionBasedStats &&
                tolower(comp_stat_string[1]) == 'u')
        {
            opt.SetUnifiedP(1);
        }
        opt.SetSmithWatermanMode(smith_waterman_value);
    }
}

void
CCompositionBasedStatsArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opt)
{
    if (args[kArgCompBasedStats])
    {
        auto_ptr<bool> ungapped(args.Exist(kArgUngapped)
                                ? new bool(args[kArgUngapped]) : 0);
        s_SetCompositionBasedStats(opt,
                                   args[kArgCompBasedStats].AsString(),
                                   args[kArgUseSWTraceback],
                                   ungapped.get());
    }

}

void
CGappedArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    // perform gapped search
#if 0
    arg_desc.AddOptionalKey(ARG_GAPPED, "gapped",
                            "Perform gapped alignment (default T, but "
                            "not available for tblastx)",
                            CArgDescriptions::eBoolean,
                            CArgDescriptions::fOptionalSeparator);
    arg_desc.AddAlias("-gapped", ARG_GAPPED);
#endif
    arg_desc.SetCurrentGroup("Extension options");
    arg_desc.AddFlag(kArgUngapped, "Perform ungapped alignment only?", true);
    arg_desc.SetCurrentGroup("");
}

void
CGappedArgs::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& options)
{
#if 0
    if (args[ARG_GAPPED] && options.GetProgram() != eTblastx)
    {
        options.SetGappedMode(args[ARG_GAPPED].AsBoolean());
    }
#endif
    options.SetGappedMode( !args[kArgUngapped] );
    if (options.GetGappedMode() == false)
        options.SetSumStatisticsMode(true);
}

void
CLargestIntronSizeArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("General search options");
    // largest intron length
    arg_desc.AddDefaultKey(kArgMaxIntronLength, "length",
                           "Length of the largest intron allowed in a translated "
                           "nucleotide sequence when linking multiple distinct "
                           "alignments (a negative value disables linking)",
                           CArgDescriptions::eInteger,
                           NStr::IntToString(kDfltArgMaxIntronLength));
    arg_desc.SetCurrentGroup("");
}

void
CLargestIntronSizeArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opt)
{
    if ( !args[kArgMaxIntronLength] )
    {
        return;
    }

    if (args[kArgMaxIntronLength].AsInteger() < 0)
    {
        opt.SetSumStatisticsMode(false);
    }
    else
    {
        opt.SetSumStatisticsMode();
        opt.SetLongestIntronLength(args[kArgMaxIntronLength].AsInteger());
    }
}

void
CFrameShiftArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("General search options");
    // applicable in blastx/tblastn, off by default
    arg_desc.AddOptionalKey(kArgFrameShiftPenalty, "frameshift",
                            "Frame shift penalty (for use with out-of-frame "
                            "gapped alignment in blastx or tblastn, default "
                            "ignored)",
                            CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgFrameShiftPenalty,
                           new CArgAllowValuesGreaterThanOrEqual(1));
    arg_desc.SetCurrentGroup("");
}

void
CFrameShiftArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opt)
{
    if (args[kArgFrameShiftPenalty])
    {
        opt.SetOutOfFrameMode();
        opt.SetFrameShiftPenalty(args[kArgFrameShiftPenalty].AsInteger());
    }
}

void
CGeneticCodeArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    if (m_Target == eQuery)
    {
        arg_desc.SetCurrentGroup("Input query options");
        // query genetic code
        arg_desc.AddDefaultKey(kArgQueryGeneticCode, "int_value",
                               "Genetic code to use to translate query",
                               CArgDescriptions::eInteger,
                               NStr::IntToString(BLAST_GENETIC_CODE));
    }
    else
    {
        arg_desc.SetCurrentGroup("General search options");
        // DB genetic code
        arg_desc.AddDefaultKey(kArgDbGeneticCode, "int_value",
                               "Genetic code to use to translate "
                               "database/subjects",
                               CArgDescriptions::eInteger,
                               NStr::IntToString(BLAST_GENETIC_CODE));
    }
    arg_desc.SetCurrentGroup("");
}

void
CGeneticCodeArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opt)
{
    const EProgram program = opt.GetProgram();

    if (m_Target == eQuery && args[kArgQueryGeneticCode])
    {
        opt.SetQueryGeneticCode(args[kArgQueryGeneticCode].AsInteger());
    }

    if (m_Target == eDatabase && args[kArgDbGeneticCode] &&
            (program == eTblastn || program == eTblastx) )
    {
        opt.SetDbGeneticCode(args[kArgDbGeneticCode].AsInteger());
    }
}

void
CGapTriggerArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("Extension options");

    const double default_value = m_QueryIsProtein
                                 ? BLAST_GAP_TRIGGER_PROT : BLAST_GAP_TRIGGER_NUCL;
    arg_desc.AddDefaultKey(kArgGapTrigger, "float_value",
                           "Number of bits to trigger gapping",
                           CArgDescriptions::eDouble,
                           NStr::DoubleToString(default_value));
    arg_desc.SetCurrentGroup("");
}

void
CGapTriggerArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opt)
{
    if (args[kArgGapTrigger])
    {
        opt.SetGapTrigger(args[kArgGapTrigger].AsDouble());
    }
}

void
CPssmEngineArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("PSSM engine options");

    // Pseudo count
    arg_desc.AddDefaultKey(kArgPSIPseudocount, "pseudocount",
                           "Pseudo-count value used when constructing PSSM",
                           CArgDescriptions::eInteger,
                           NStr::IntToString(PSI_PSEUDO_COUNT_CONST));

    if (m_IsDeltaBlast)
    {
        arg_desc.AddDefaultKey(kArgDomainInclusionEThreshold, "ethresh",
                               "E-value inclusion threshold for alignments "
                               "with conserved domains",
                               CArgDescriptions::eDouble,
                               NStr::DoubleToString(DELTA_INCLUSION_ETHRESH));
    }

    // Evalue inclusion threshold
    arg_desc.AddDefaultKey(kArgPSIInclusionEThreshold, "ethresh",
                           "E-value inclusion threshold for pairwise alignments",
                           CArgDescriptions::eDouble,
                           NStr::DoubleToString(PSI_INCLUSION_ETHRESH));

    arg_desc.SetCurrentGroup("");
}

void
CPssmEngineArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opt)
{
    if (args[kArgPSIPseudocount])
    {
        opt.SetPseudoCount(args[kArgPSIPseudocount].AsInteger());
    }

    if (args[kArgPSIInclusionEThreshold])
    {
        opt.SetInclusionThreshold(args[kArgPSIInclusionEThreshold].AsDouble());
    }

    if (args.Exist(kArgDomainInclusionEThreshold)
            && args[kArgDomainInclusionEThreshold])
    {

        opt.SetDomainInclusionThreshold(
            args[kArgDomainInclusionEThreshold].AsDouble());
    }
}

void
CPsiBlastArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{

    if (m_DbTarget == eNucleotideDb)
    {
        arg_desc.SetCurrentGroup("PSI-TBLASTN options");

        // PSI-tblastn checkpoint
        arg_desc.AddOptionalKey(kArgPSIInputChkPntFile, "psi_chkpt_file",
                                "PSI-TBLASTN checkpoint file",
                                CArgDescriptions::eInputFile);
        arg_desc.SetDependency(kArgPSIInputChkPntFile,
                               CArgDescriptions::eExcludes,
                               kArgRemote);
    }
    else
    {
        arg_desc.SetCurrentGroup("PSI-BLAST options");

        // Number of iterations
        arg_desc.AddDefaultKey(kArgPSINumIterations, "int_value",
                               "Number of iterations to perform",
                               CArgDescriptions::eInteger,
                               NStr::IntToString(1));
        arg_desc.SetConstraint(kArgPSINumIterations,
                               new CArgAllowValuesGreaterThanOrEqual(1));
        arg_desc.SetDependency(kArgPSINumIterations,
                               CArgDescriptions::eExcludes,
                               kArgRemote);
        // checkpoint file
        arg_desc.AddOptionalKey(kArgPSIOutputChkPntFile, "checkpoint_file",

                                "File name to store checkpoint file",
                                CArgDescriptions::eOutputFile);
        // ASCII matrix file
        arg_desc.AddOptionalKey(kArgAsciiPssmOutputFile, "ascii_mtx_file",
                                "File name to store ASCII version of PSSM",
                                CArgDescriptions::eOutputFile);

        if (!m_IsDeltaBlast)
        {
            vector<string> msa_exclusions;
            msa_exclusions.push_back(kArgPSIInputChkPntFile);
            msa_exclusions.push_back(kArgQuery);
            msa_exclusions.push_back(kArgQueryLocation);
            // pattern and MSA is not supported
            msa_exclusions.push_back(kArgPHIPatternFile);
            arg_desc.SetCurrentGroup("");
            arg_desc.SetCurrentGroup("");

            // MSA restart file
            arg_desc.SetCurrentGroup("PSSM engine options");
            arg_desc.AddOptionalKey(kArgMSAInputFile, "align_restart",
                                    "File name of multiple sequence alignment to "
                                    "restart PSI-BLAST",
                                    CArgDescriptions::eInputFile);
            ITERATE(vector<string>, exclusion, msa_exclusions)
            {
                arg_desc.SetDependency(kArgMSAInputFile,
                                       CArgDescriptions::eExcludes,
                                       *exclusion);
            }

            arg_desc.AddOptionalKey(kArgMSAMasterIndex, "index",
                                    "Ordinal number (1-based index) of the sequence"
                                    " to use as a master in the multiple sequence "
                                    "alignment. If not provided, the first sequence"
                                    " in the multiple sequence alignment will be "
                                    "used", CArgDescriptions::eInteger);
            arg_desc.SetConstraint(kArgMSAMasterIndex,
                                   new CArgAllowValuesGreaterThanOrEqual(1));
            ITERATE(vector<string>, exclusion, msa_exclusions)
            {
                arg_desc.SetDependency(kArgMSAMasterIndex,
                                       CArgDescriptions::eExcludes,
                                       *exclusion);
            }
            arg_desc.SetDependency(kArgMSAMasterIndex,
                                   CArgDescriptions::eRequires,
                                   kArgMSAInputFile);
            arg_desc.SetDependency(kArgMSAMasterIndex,
                                   CArgDescriptions::eExcludes,
                                   kArgIgnoreMsaMaster);

            arg_desc.AddFlag(kArgIgnoreMsaMaster,
                             "Ignore the master sequence when creating PSSM", true);

            vector<string> ignore_pssm_master_exclusions;
            ignore_pssm_master_exclusions.push_back(kArgMSAMasterIndex);
            ignore_pssm_master_exclusions.push_back(kArgPSIInputChkPntFile);
            ignore_pssm_master_exclusions.push_back(kArgQuery);
            ignore_pssm_master_exclusions.push_back(kArgQueryLocation);
            ITERATE(vector<string>, exclusion, msa_exclusions)
            {
                arg_desc.SetDependency(kArgIgnoreMsaMaster,
                                       CArgDescriptions::eExcludes,
                                       *exclusion);
            }
            arg_desc.SetDependency(kArgIgnoreMsaMaster,
                                   CArgDescriptions::eRequires,
                                   kArgMSAInputFile);

            // PSI-BLAST checkpoint
            arg_desc.AddOptionalKey(kArgPSIInputChkPntFile, "psi_chkpt_file",
                                    "PSI-BLAST checkpoint file",
                                    CArgDescriptions::eInputFile);
        }
    }

    if (!m_IsDeltaBlast)
    {
        arg_desc.SetDependency(kArgPSIInputChkPntFile,
                               CArgDescriptions::eExcludes,
                               kArgQuery);
        arg_desc.SetDependency(kArgPSIInputChkPntFile,
                               CArgDescriptions::eExcludes,
                               kArgQueryLocation);
    }
    arg_desc.SetCurrentGroup("");
}

CRef<CPssmWithParameters>
CPsiBlastArgs::x_CreatePssmFromMsa(CNcbiIstream& input_stream,
                                   CBlastOptions& opt, bool save_ascii_pssm,
                                   unsigned int msa_master_idx,
                                   bool ignore_pssm_tmplt_seq)
{
    // FIXME get these from CBlastOptions
    CPSIBlastOptions psiblast_opts;
    PSIBlastOptionsNew(&psiblast_opts);
    psiblast_opts->nsg_compatibility_mode = ignore_pssm_tmplt_seq;

    CPSIDiagnosticsRequest diags(PSIDiagnosticsRequestNewEx(save_ascii_pssm));
    CPsiBlastInputClustalW pssm_input(input_stream, *psiblast_opts,
                                      opt.GetMatrixName(), diags, NULL, 0,
                                      opt.GetGapOpeningCost(),
                                      opt.GetGapExtensionCost(),
                                      msa_master_idx);
    CPssmEngine pssm_engine(&pssm_input);
    return pssm_engine.Run();
}

void
CPsiBlastArgs::ExtractAlgorithmOptions(const CArgs& args,
                                       CBlastOptions& opt)
{
    if (m_DbTarget == eProteinDb)
    {
        if (args[kArgPSINumIterations])
        {
            if(m_NumIterations == 1)
                m_NumIterations = args[kArgPSINumIterations].AsInteger();
        }
        if (args.Exist(kArgPSIOutputChkPntFile) &&
                args[kArgPSIOutputChkPntFile])
        {
            m_CheckPointOutput.Reset
            (new CAutoOutputFileReset
             (args[kArgPSIOutputChkPntFile].AsString()));
        }
        const bool kSaveAsciiPssm = args[kArgAsciiPssmOutputFile];
        if (kSaveAsciiPssm)
        {
            m_AsciiMatrixOutput.Reset
            (new CAutoOutputFileReset
             (args[kArgAsciiPssmOutputFile].AsString()));
        }
        if (args.Exist(kArgMSAInputFile) && args[kArgMSAInputFile])
        {
            CNcbiIstream& in = args[kArgMSAInputFile].AsInputFile();
            unsigned int msa_master_idx = 0;
            if (args[kArgMSAMasterIndex])
            {
                msa_master_idx = args[kArgMSAMasterIndex].AsInteger() - 1;
            }
            m_Pssm = x_CreatePssmFromMsa(in, opt, kSaveAsciiPssm,
                                         msa_master_idx,
                                         args[kArgIgnoreMsaMaster]);
        }
        if (!m_IsDeltaBlast)
        {
            opt.SetIgnoreMsaMaster(args[kArgIgnoreMsaMaster]);
        }
    }

    if (args.Exist(kArgPSIInputChkPntFile) && args[kArgPSIInputChkPntFile])
    {
        CNcbiIstream& in = args[kArgPSIInputChkPntFile].AsInputFile();
        _ASSERT(m_Pssm.Empty());
        m_Pssm.Reset(new CPssmWithParameters);
        try
        {
            switch (CFormatGuess().Format(in))
            {
            case CFormatGuess::eBinaryASN:
                in >> MSerial_AsnBinary >> *m_Pssm;
                break;
            case CFormatGuess::eTextASN:
                in >> MSerial_AsnText >> *m_Pssm;
                break;
            case CFormatGuess::eXml:
                in >> MSerial_Xml >> *m_Pssm;
                break;
            default:
                NCBI_THROW(CInputException, eInvalidInput,
                           "Unsupported format for PSSM");
            }
        }
        catch (const CSerialException&)
        {
            string msg("Unrecognized format for PSSM in ");
            msg += args[kArgPSIInputChkPntFile].AsString() + " (must be ";
            msg += "PssmWithParameters)";
            NCBI_THROW(CInputException, eInvalidInput, msg);
        }
        _ASSERT(m_Pssm.NotEmpty());
    }
}

void
CPhiBlastArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("PHI-BLAST options");

    arg_desc.AddOptionalKey(kArgPHIPatternFile, "file",
                            "File name containing pattern to search",
                            CArgDescriptions::eInputFile);
    arg_desc.SetDependency(kArgPHIPatternFile,
                           CArgDescriptions::eExcludes,
                           kArgPSIInputChkPntFile);

    arg_desc.SetCurrentGroup("");
}

void
CPhiBlastArgs::ExtractAlgorithmOptions(const CArgs& args,
                                       CBlastOptions& opt)
{
    if (args.Exist(kArgPHIPatternFile) && args[kArgPHIPatternFile])
    {
        CNcbiIstream& in = args[kArgPHIPatternFile].AsInputFile();
        in.clear();
        in.seekg(0);
        char buffer[4096];
        string line;
        string pattern;
        string name;
        while (in.getline(buffer, 4096))
        {
            line = buffer;
            string ltype = line.substr(0, 2);
            if (ltype == "ID")
                name = line.substr(5);
            else if (ltype == "PA")
                pattern = line.substr(5);
        }
        if (!pattern.empty())
            opt.SetPHIPattern(pattern.c_str(),
                              (Blast_QueryIsNucleotide(opt.GetProgramType())
                               ? true : false));
        else
            NCBI_THROW(CInputException, eInvalidInput,
                       "PHI pattern not read");
    }
}

void
CDeltaBlastArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("DELTA-BLAST options");

    arg_desc.AddDefaultKey(kArgRpsDb, "database_name", "BLAST domain "
                           "database name", CArgDescriptions::eString,
                           kDfltArgRpsDb);

    arg_desc.AddFlag(kArgShowDomainHits, "Show domain hits");
    arg_desc.SetDependency(kArgShowDomainHits, CArgDescriptions::eExcludes,
                           kArgRemote);
}

void
CDeltaBlastArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opt)
{
    m_DomainDb.Reset(new CSearchDatabase(args[kArgRpsDb].AsString(),
                                         CSearchDatabase::eBlastDbIsProtein));

    if (args.Exist(kArgShowDomainHits))
    {
        m_ShowDomainHits = args[kArgShowDomainHits];
    }
}

void
CIgBlastArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("Ig-BLAST options");
    const static string suffix = "VDJ";
    const static int df_num_align[3] = {3,3,3};
    int num_genes = (m_IsProtein) ? 1 : 3;

    for (int gene=0; gene<num_genes; ++gene)
    {
        // Subject sequence input
        /* TODO disabled for now
        string arg_sub = kArgGLSubject;
        arg_sub.push_back(suffix[gene]);
        arg_desc.AddOptionalKey(arg_sub , "filename",
                            "Germline subject sequence to align",
                            CArgDescriptions::eInputFile);
        */
        // Germline database file name
        string arg_db = kArgGLDatabase;
        arg_db.push_back(suffix[gene]);
        arg_desc.AddOptionalKey(arg_db, "germline_database_name",
                                "Germline database name",
                                CArgDescriptions::eString);
        //arg_desc.SetDependency(arg_db, CArgDescriptions::eExcludes, arg_sub);
        // Number of alignments to show
        string arg_na = kArgGLNumAlign;
        arg_na.push_back(suffix[gene]);
        arg_desc.AddDefaultKey(arg_na, "int_value",
                               "Number of Germline sequences to show alignments for",
                               CArgDescriptions::eInteger,
                               NStr::IntToString(df_num_align[gene]));
        //arg_desc.SetConstraint(arg_na,
        //                   new CArgAllowValuesBetween(0, 4));
        // Seqidlist
        arg_desc.AddOptionalKey(arg_db + "_seqidlist", "filename",
                                "Restrict search of germline database to list of SeqIds's",
                                CArgDescriptions::eString);
    }

    if (!m_IsProtein)
    {
        arg_desc.AddOptionalKey(kArgGLChainType, "filename",
                                "File containing the coding frame start positions for sequences in germline J database",
                                CArgDescriptions::eString);

        arg_desc.AddOptionalKey(kArgMinDMatch, "min_D_match",
                                "Required minimal number of D gene matches ",
                                CArgDescriptions::eInteger);
        arg_desc.SetConstraint(kArgMinDMatch,
                               new CArgAllowValuesGreaterThanOrEqual(5));

    }

    arg_desc.AddDefaultKey(kArgGLOrigin, "germline_origin",
                           "The organism for your query sequence (i.e., human, mouse, etc.)",
                           CArgDescriptions::eString, "human");
    //arg_desc.SetConstraint(kArgGLOrigin, &(*new CArgAllow_Strings, "human", "mouse"));

    arg_desc.AddDefaultKey(kArgGLDomainSystem, "domain_system",
                           "Domain system to be used for segment annotation",
                           CArgDescriptions::eString, "kabat");
    arg_desc.SetConstraint(kArgGLDomainSystem, &(*new CArgAllow_Strings, "kabat", "imgt"));

    arg_desc.AddFlag(kArgGLFocusV, "Should the search only be for V segment (effective only for non-germline database search using -db option)?", true);

    if (! m_IsProtein)
    {
        arg_desc.AddFlag(kArgTranslate, "Show translated alignments", true);
    }

    arg_desc.SetCurrentGroup("");
}

static string s_RegisterOMDataLoader(CRef<CSeqDB> db_handle)
{
    // the blast formatter requires that the database coexist in
    // the same scope with the query sequences
    CRef<CObjectManager> om = CObjectManager::GetInstance();
    CBlastDbDataLoader::RegisterInObjectManager(*om, db_handle, true,
            CObjectManager::eDefault,
            CBlastDatabaseArgs::kSubjectsDataLoaderPriority);
    CBlastDbDataLoader::SBlastDbParam param(db_handle);
    string retval(CBlastDbDataLoader::GetLoaderNameFromArgs(param));
    _TRACE("Registering " << retval << " at priority " <<
           (int)CBlastDatabaseArgs::kSubjectsDataLoaderPriority);
    return retval;
}

void
CIgBlastArgs::ExtractAlgorithmOptions(const CArgs& args,
                                      CBlastOptions& opts)
{
    string paths[3];
    CNcbiEnvironment env;
    CMetaRegistry::SEntry sentry =
        CMetaRegistry::Load("ncbi", CMetaRegistry::eName_RcOrIni);
    paths[0] = CDirEntry::NormalizePath(CDir::GetCwd(), eFollowLinks);
    paths[1] = CDirEntry::NormalizePath(env.Get("IGDATA"), eFollowLinks);
    if (sentry.registry)
    {
        paths[2] = CDirEntry::NormalizePath(sentry.registry->Get("BLAST","IGDATA"), eFollowLinks);
    }
    else
    {
#if defined(NCBI_OS_DARWIN)
        paths[2] = "/usr/local/ncbi/igblast/data";
#else
        paths[2] = paths[0];
#endif
    }

    m_IgOptions.Reset(new CIgBlastOptions());

    CBlastDatabaseArgs::EMoleculeType mol_type = Blast_SubjectIsNucleotide(opts.GetProgramType())
            ? CSearchDatabase::eBlastDbIsNucleotide
            : CSearchDatabase::eBlastDbIsProtein;

    m_IgOptions->m_IsProtein = (mol_type == CSearchDatabase::eBlastDbIsProtein);
    m_IgOptions->m_Origin = args[kArgGLOrigin].AsString();
    m_IgOptions->m_DomainSystem = args[kArgGLDomainSystem].AsString();
    m_IgOptions->m_FocusV = args.Exist(kArgGLFocusV) ? args[kArgGLFocusV] : false;
    m_IgOptions->m_Translate = args.Exist(kArgTranslate) ? args[kArgTranslate] : false;
    if (!m_IsProtein)
    {
        string aux_file = (args.Exist(kArgGLChainType) && args[kArgGLChainType])
                          ? args[kArgGLChainType].AsString()
                          : m_IgOptions->m_Origin + "_gl.aux";
        m_IgOptions->m_AuxFilename = aux_file;
        for (int i=0; i<3; i++)
        {
            string aux_path = CDirEntry::ConcatPath(paths[i], aux_file);
            CDirEntry entry(aux_path);
            if (entry.Exists() && entry.IsFile())
            {
                m_IgOptions->m_AuxFilename = aux_path;
                break;
            }
        }
    }

    _ASSERT(m_IsProtein == m_IgOptions->m_IsProtein);

    m_Scope.Reset(new CScope(*CObjectManager::GetInstance()));

    // default germline database name for annotation
    for (int i=0; i<3; i++)
    {
        string int_data = CDirEntry::ConcatPath(paths[i], "internal_data");
        CDirEntry entry(int_data);
        if (entry.Exists() && entry.IsDir())
        {
            m_IgOptions->m_IgDataPath = int_data;
            break;
        }
    }

    string df_db_name = CDirEntry::ConcatPath(
                            CDirEntry::ConcatPath(m_IgOptions->m_IgDataPath,
                                    m_IgOptions->m_Origin), m_IgOptions->m_Origin + "_V");
    CRef<CSearchDatabase> db(new CSearchDatabase(df_db_name, mol_type));
    m_IgOptions->m_Db[3].Reset(new CLocalDbAdapter(*db));
    try
    {
        db->GetSeqDb();
    }
    catch(...)
    {
        NCBI_THROW(CInputException, eInvalidInput,
                   "Germline annotation database " + df_db_name + " could not be found in [internal_data] directory");
    }

    m_IgOptions->m_Min_D_match = 5;
    if (args.Exist(kArgMinDMatch) && args[kArgMinDMatch])
    {
        m_IgOptions->m_Min_D_match = args[kArgMinDMatch].AsInteger();
    }


    CRef<CBlastOptionsHandle> opts_hndl;
    if (m_IgOptions->m_IsProtein)
    {
        opts_hndl.Reset(CBlastOptionsFactory::Create(eBlastp));
    }
    else
    {
        opts_hndl.Reset(CBlastOptionsFactory::Create(eBlastn));
    }

    const static string suffix = "VDJ";
    int num_genes = (m_IsProtein) ? 1: 3;
    for (int gene=0; gene< num_genes; ++gene)
    {
        string arg_sub = kArgGLSubject;
        string arg_db  = kArgGLDatabase;
        string arg_na  = kArgGLNumAlign;

        arg_sub.push_back(suffix[gene]);
        arg_db.push_back(suffix[gene]);
        arg_na.push_back(suffix[gene]);

        m_IgOptions->m_NumAlign[gene] = args[arg_na].AsInteger();

        if (args.Exist(arg_sub) && args[arg_sub])
        {
            CNcbiIstream& subj_input_stream = args[arg_sub].AsInputFile();
            TSeqRange subj_range;

            const bool parse_deflines = args.Exist(kArgParseDeflines)
                                        ? bool(args[kArgParseDeflines])
                                        : kDfltArgParseDeflines;
            const bool use_lcase_masks = args.Exist(kArgUseLCaseMasking)
                                         ? bool(args[kArgUseLCaseMasking])
                                         : kDfltArgUseLCaseMasking;
            CRef<blast::CBlastQueryVector> subjects;
            CRef<CScope> scope = ReadSequencesToBlast(subj_input_stream,
                                 m_IgOptions->m_IsProtein,
                                 subj_range, parse_deflines,
                                 use_lcase_masks, subjects);
            m_Scope->AddScope(*scope,
                              CBlastDatabaseArgs::kSubjectsDataLoaderPriority);
            CRef<IQueryFactory> sub_seqs(
                new blast::CObjMgr_QueryFactory(*subjects));
            m_IgOptions->m_Db[gene].Reset(new CLocalDbAdapter(
                                              sub_seqs, opts_hndl));
        }
        else
        {
            string gl_db_name = m_IgOptions->m_Origin + "_gl_";
            gl_db_name.push_back(suffix[gene]);
            string db_name = (args.Exist(arg_db) && args[arg_db])
                             ?  args[arg_db].AsString() :  gl_db_name;
            db.Reset(new CSearchDatabase(db_name, mol_type));

            if (args.Exist(arg_db + "_seqidlist") && args[arg_db + "_seqidlist"])
            {
                string fn(SeqDB_ResolveDbPath(args[arg_db + "_seqidlist"].AsString()));
                db->SetGiList(CRef<CSeqDBGiList> (new CSeqDBFileGiList(fn,
                                                  CSeqDBFileGiList::eSiList)));
            }

            m_IgOptions->m_Db[gene].Reset(new CLocalDbAdapter(*db));
            m_Scope->AddDataLoader(s_RegisterOMDataLoader(db->GetSeqDb()));
        }
    }
}

void
CQueryOptionsArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{

    arg_desc.SetCurrentGroup("Query filtering options");
    // lowercase masking
    arg_desc.AddFlag(kArgUseLCaseMasking,
                     "Use lower case filtering in query and subject sequence(s)?", true);

    arg_desc.SetCurrentGroup("Input query options");
    // query location
    arg_desc.AddOptionalKey(kArgQueryLocation, "range",
                            "Location on the query sequence in 1-based offsets "
                            "(Format: start-stop)",
                            CArgDescriptions::eString);

    if ( !m_QueryCannotBeNucl )
    {
        // search strands
        arg_desc.AddDefaultKey(kArgStrand, "strand",
                               "Query strand(s) to search against database/subject",
                               CArgDescriptions::eString, kDfltArgStrand);
        arg_desc.SetConstraint(kArgStrand, &(*new CArgAllow_Strings,
                                             kDfltArgStrand, "plus", "minus"));
    }

    arg_desc.SetCurrentGroup("Miscellaneous options");
    arg_desc.AddFlag(kArgParseDeflines,
                     "Should the query and subject defline(s) be parsed?", true);

    arg_desc.SetCurrentGroup("");
}

void
CQueryOptionsArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opt)
{
    // Get the strand
    {
        m_Strand = eNa_strand_unknown;

        if (!Blast_QueryIsProtein(opt.GetProgramType()) && args[kArgStrand])
        {
            const string& kStrand = args[kArgStrand].AsString();
            if (kStrand == "both")
            {
                m_Strand = eNa_strand_both;
            }
            else if (kStrand == "plus")
            {
                m_Strand = eNa_strand_plus;
            }
            else if (kStrand == "minus")
            {
                m_Strand = eNa_strand_minus;
            }
            else
            {
                abort();
            }
        }
    }

    // set the sequence range
    if (args[kArgQueryLocation])
    {
        m_Range = ParseSequenceRange(args[kArgQueryLocation].AsString(),
                                     "Invalid specification of query location");
    }

    m_UseLCaseMask = static_cast<bool>(args[kArgUseLCaseMasking]);
    m_ParseDeflines = static_cast<bool>(args[kArgParseDeflines]);
}

CBlastDatabaseArgs::CBlastDatabaseArgs(bool request_mol_type /* = false */,
                                       bool is_rpsblast /* = false */,
                                       bool is_igblast  /* = false */,
                                       bool is_deltablast /* = false */)
    : m_RequestMoleculeType(request_mol_type),
    m_IsRpsBlast(is_rpsblast),
    m_IsIgBlast(is_igblast),
    m_IsDeltaBlast(is_deltablast),
    m_IsProtein(true),
    m_SupportsDatabaseMasking(false)
{}

bool
CBlastDatabaseArgs::HasBeenSet(const CArgs& args)
{
    if ( (args.Exist(kArgDb) && args[kArgDb].HasValue()) ||
            (args.Exist(kArgSubject) && args[kArgSubject].HasValue()) )
    {
        return true;
    }
    return false;
}

void
CBlastDatabaseArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("General search options");
    // database filename
    arg_desc.AddOptionalKey(kArgDb, "database_name", "BLAST database name",
                            CArgDescriptions::eString);
    arg_desc.SetCurrentGroup("");

    if (m_RequestMoleculeType)
    {
        arg_desc.AddKey(kArgDbType, "database_type",
                        "BLAST database molecule type",
                        CArgDescriptions::eString);
        arg_desc.SetConstraint(kArgDbType,
                               &(*new CArgAllow_Strings, "prot", "nucl"));
    }

    vector<string> database_args;
    database_args.push_back(kArgDb);
    database_args.push_back(kArgGiList);
    database_args.push_back(kArgSeqIdList);
    database_args.push_back(kArgNegativeGiList);
    if (m_SupportsDatabaseMasking)
    {
        database_args.push_back(kArgDbSoftMask);
        database_args.push_back(kArgDbHardMask);
    }

    // DB size
    arg_desc.SetCurrentGroup("Statistical options");
    arg_desc.AddOptionalKey(kArgDbSize, "num_letters",
                            "Effective length of the database ",
                            CArgDescriptions::eInt8);

    arg_desc.SetCurrentGroup("Restrict search or results");
    // GI list
    arg_desc.AddOptionalKey(kArgGiList, "filename",
                            "Restrict search of database to list of GI's",
                            CArgDescriptions::eString);
    // SeqId list
    arg_desc.AddOptionalKey(kArgSeqIdList, "filename",
                            "Restrict search of database to list of SeqId's",
                            CArgDescriptions::eString);
    // Negative GI list
    arg_desc.AddOptionalKey(kArgNegativeGiList, "filename",
                            "Restrict search of database to everything except the listed GIs",
                            CArgDescriptions::eString);
    arg_desc.SetDependency(kArgGiList, CArgDescriptions::eExcludes,
                           kArgNegativeGiList);
    arg_desc.SetDependency(kArgGiList, CArgDescriptions::eExcludes,
                           kArgSeqIdList);
    arg_desc.SetDependency(kArgSeqIdList, CArgDescriptions::eExcludes,
                           kArgNegativeGiList);
    // Entrez Query
    // Deltablast does not support remote blast
    if (!m_IsDeltaBlast)
    {
        arg_desc.AddOptionalKey(kArgEntrezQuery, "entrez_query",
                                "Restrict search with the given Entrez query",
                                CArgDescriptions::eString);


        // For now, disable pairing -remote with either -gilist or
        // -negative_gilist as this is not implemented in the BLAST server
        arg_desc.SetDependency(kArgGiList, CArgDescriptions::eExcludes,
                               kArgRemote);
        arg_desc.SetDependency(kArgSeqIdList, CArgDescriptions::eExcludes,
                               kArgRemote);
        arg_desc.SetDependency(kArgNegativeGiList, CArgDescriptions::eExcludes,
                               kArgRemote);

        // Entrez query currently requires the -remote option
        arg_desc.SetDependency(kArgEntrezQuery, CArgDescriptions::eRequires,
                               kArgRemote);
    }

#if ((!defined(NCBI_COMPILER_WORKSHOP) || (NCBI_COMPILER_VERSION  > 550)) && \
     (!defined(NCBI_COMPILER_MIPSPRO)) )
    // Masking of database
    if (m_SupportsDatabaseMasking)
    {
        arg_desc.AddOptionalKey(kArgDbSoftMask,
                                "filtering_algorithm",
                                "Filtering algorithm ID to apply to the BLAST database as soft "
                                "masking",
                                CArgDescriptions::eString);
        arg_desc.SetDependency(kArgDbSoftMask, CArgDescriptions::eExcludes,
                               kArgDbHardMask);

        arg_desc.AddOptionalKey(kArgDbHardMask,
                                "filtering_algorithm",
                                "Filtering algorithm ID to apply to the BLAST database as hard "
                                "masking",
                                CArgDescriptions::eString);
    }
#endif

    // There is no RPS-BLAST 2 sequences
    if ( !m_IsRpsBlast )
    {
        arg_desc.SetCurrentGroup("BLAST-2-Sequences options");
        // subject sequence input (for bl2seq)
        arg_desc.AddOptionalKey(kArgSubject, "subject_input_file",
                                "Subject sequence(s) to search",
                                CArgDescriptions::eInputFile);
        ITERATE(vector<string>, dbarg, database_args)
        {
            arg_desc.SetDependency(kArgSubject, CArgDescriptions::eExcludes,
                                   *dbarg);
        }

        // subject location
        arg_desc.AddOptionalKey(kArgSubjectLocation, "range",
                                "Location on the subject sequence in 1-based offsets "
                                "(Format: start-stop)",
                                CArgDescriptions::eString);
        ITERATE(vector<string>, dbarg, database_args)
        {
            arg_desc.SetDependency(kArgSubjectLocation,
                                   CArgDescriptions::eExcludes,
                                   *dbarg);
        }
        // Because Blast4-subject does not support Seq-locs, specifying a
        // subject range does not work for remote searches
        arg_desc.SetDependency(kArgSubjectLocation,
                               CArgDescriptions::eExcludes, kArgRemote);
    }

    arg_desc.SetCurrentGroup("");
}

void
CBlastDatabaseArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opts)
{
    EMoleculeType mol_type = Blast_SubjectIsNucleotide(opts.GetProgramType())
                             ? CSearchDatabase::eBlastDbIsNucleotide
                             : CSearchDatabase::eBlastDbIsProtein;
    m_IsProtein = (mol_type == CSearchDatabase::eBlastDbIsProtein);

    if (args.Exist(kArgDb) && args[kArgDb])
    {

        m_SearchDb.Reset(new CSearchDatabase(args[kArgDb].AsString(),
                                             mol_type));

        if (args.Exist(kArgGiList) && args[kArgGiList])
        {
            string fn(SeqDB_ResolveDbPath(args[kArgGiList].AsString()));
            m_SearchDb->SetGiList(CRef<CSeqDBGiList> (new CSeqDBFileGiList(fn)));

        }
        else if (args.Exist(kArgNegativeGiList) && args[kArgNegativeGiList])
        {
            string fn(SeqDB_ResolveDbPath(args[kArgNegativeGiList].AsString()));
            m_SearchDb->SetNegativeGiList(CRef<CSeqDBGiList> (new CSeqDBFileGiList(fn)));

        }
        else if (args.Exist(kArgSeqIdList) && args[kArgSeqIdList])
        {
            string fn(SeqDB_ResolveDbPath(args[kArgSeqIdList].AsString()));
            m_SearchDb->SetGiList(CRef<CSeqDBGiList> (new CSeqDBFileGiList(fn,
                                  CSeqDBFileGiList::eSiList)));
        }

        if (args.Exist(kArgEntrezQuery) && args[kArgEntrezQuery])
            m_SearchDb->SetEntrezQueryLimitation(args[kArgEntrezQuery].AsString());

#if ((!defined(NCBI_COMPILER_WORKSHOP) || (NCBI_COMPILER_VERSION  > 550)) && \
     (!defined(NCBI_COMPILER_MIPSPRO)) )
        if (args.Exist(kArgDbSoftMask) && args[kArgDbSoftMask])
        {
            m_SearchDb->SetFilteringAlgorithm(args[kArgDbSoftMask].AsString(), eSoftSubjMasking);
        }
        else if (args.Exist(kArgDbHardMask) && args[kArgDbHardMask])
        {
            m_SearchDb->SetFilteringAlgorithm(args[kArgDbHardMask].AsString(), eHardSubjMasking);
        }
#endif
    }
    else if (args.Exist(kArgSubject) && args[kArgSubject])
    {

        CNcbiIstream& subj_input_stream = args[kArgSubject].AsInputFile();
        TSeqRange subj_range;
        if (args.Exist(kArgSubjectLocation) && args[kArgSubjectLocation])
        {
            subj_range =
                ParseSequenceRange(args[kArgSubjectLocation].AsString(),
                                   "Invalid specification of subject location");
        }

        const bool parse_deflines = args.Exist(kArgParseDeflines)
                                    ? bool(args[kArgParseDeflines])
                                    : kDfltArgParseDeflines;
        const bool use_lcase_masks = args.Exist(kArgUseLCaseMasking)
                                     ? bool(args[kArgUseLCaseMasking])
                                     : kDfltArgUseLCaseMasking;
        CRef<blast::CBlastQueryVector> subjects;
        m_Scope = ReadSequencesToBlast(subj_input_stream, IsProtein(),
                                       subj_range, parse_deflines,
                                       use_lcase_masks, subjects);
        m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects));

    }
    else if (!m_IsIgBlast)
    {
        // IgBlast permits use of germline database
        NCBI_THROW(CInputException, eInvalidInput,
                   "Either a BLAST database or subject sequence(s) must be specified");
    }

    if (opts.GetEffectiveSearchSpace() != 0)
    {
        // no need to set any other options, as this trumps them
        return;
    }

    if (args[kArgDbSize])
    {
        opts.SetDbLength(args[kArgDbSize].AsInt8());
    }

}

void
CFormattingArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("Formatting options");

    const string kOutputFormatDescription = string(
            "alignment view options:\n"
            "  0 = pairwise,\n"
            "  1 = query-anchored showing identities,\n"
            "  2 = query-anchored no identities,\n"
            "  3 = flat query-anchored, show identities,\n"
            "  4 = flat query-anchored, no identities,\n"
            "  5 = XML Blast output,\n"
            "  6 = tabular,\n"
            "  7 = tabular with comment lines,\n"
            "  8 = Text ASN.1,\n"
            "  9 = Binary ASN.1,\n"
            " 10 = Comma-separated values,\n"
            " 11 = BLAST archive format (ASN.1) \n\n"
            "Options 6, 7, and 10 can be additionally configured to produce\n"
            "a custom format specified by space delimited format specifiers.\n"
            "The supported format specifiers are:\n") +
                                            DescribeTabularOutputFormatSpecifiers() +
                                            string("\n");

    const string kIgOutputFormatDescription = string(
                "alignment view options:\n"
                "  3 = flat query-anchored, show identities,\n"
                "  4 = flat query-anchored, no identities,\n"
                "  7 = tabular with comment lines\n\n"
                "Options 7 can be additionally configured to produce\n"
                "a custom format specified by space delimited format specifiers.\n"
                "The supported format specifiers are:\n") +
            DescribeTabularOutputFormatSpecifiers() +
            string("\n");

    int dft_outfmt = (m_IsIgBlast) ? 3 : kDfltArgOutputFormat;

    // alignment view
    arg_desc.AddDefaultKey(kArgOutputFormat, "format",
                           (m_IsIgBlast) ? kIgOutputFormatDescription
                           : kOutputFormatDescription,
                           CArgDescriptions::eString,
                           NStr::IntToString(dft_outfmt));

    // show GIs in deflines
    arg_desc.AddFlag(kArgShowGIs, "Show NCBI GIs in deflines?", true);

    // number of one-line descriptions to display
    arg_desc.AddOptionalKey(kArgNumDescriptions, "int_value",
                            "Number of database sequences to show one-line "
                            "descriptions for\n"
                            "Not applicable for outfmt > 4\n"
                            "Default = `"+ NStr::IntToString(m_DfltNumDescriptions)+ "'",
                            CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgNumDescriptions,
                           new CArgAllowValuesGreaterThanOrEqual(0));

    // number of alignments per DB sequence
    arg_desc.AddOptionalKey(kArgNumAlignments, "int_value",
                            "Number of database sequences to show alignments for\n"
                            "Default = `" + NStr::IntToString(m_DfltNumAlignments) + "'",
                            CArgDescriptions::eInteger );
    arg_desc.SetConstraint(kArgNumAlignments,
                           new CArgAllowValuesGreaterThanOrEqual(0));

    // Produce HTML?
    if(!m_IsIgBlast)
    {
        arg_desc.AddFlag(kArgProduceHtml, "Produce HTML output?", true);
    }
    /// Hit list size, listed here for convenience only
    arg_desc.SetCurrentGroup("Restrict search or results");
    arg_desc.AddOptionalKey(kArgMaxTargetSequences, "num_sequences",
                            "Maximum number of aligned sequences to keep \n"
                            "Not applicable for outfmt <= 4\n"
                            "Default = `" + NStr::IntToString(BLAST_HITLIST_SIZE) + "'",
                            CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMaxTargetSequences,
                           new CArgAllowValuesGreaterThanOrEqual(1));
    arg_desc.SetDependency(kArgMaxTargetSequences,
                           CArgDescriptions::eExcludes,
                           kArgNumDescriptions);
    arg_desc.SetDependency(kArgMaxTargetSequences,
                           CArgDescriptions::eExcludes,
                           kArgNumAlignments);

    arg_desc.SetCurrentGroup("");
}

bool
CFormattingArgs::ArchiveFormatRequested(const CArgs& args) const
{
    EOutputFormat output_fmt;
    string ignore;
    ParseFormattingString(args, output_fmt, ignore);
    return (output_fmt == eArchiveFormat ? true : false);
}

void
CFormattingArgs::ParseFormattingString(const CArgs& args,
                                       EOutputFormat& fmt_type,
                                       string& custom_fmt_spec) const
{
    custom_fmt_spec.clear();
    if (args[kArgOutputFormat])
    {
        string fmt_choice =
            NStr::TruncateSpaces(args[kArgOutputFormat].AsString());
        string::size_type pos;
        if ( (pos = fmt_choice.find_first_of(' ')) != string::npos)
        {
            custom_fmt_spec.assign(fmt_choice, pos+1,
                                   fmt_choice.size()-(pos+1));
            fmt_choice.erase(pos);
        }
        int val = 0;
        try
        {
            val = NStr::StringToInt(fmt_choice);
        }
        catch (const CStringException&)     // probably a conversion error
        {
            CNcbiOstrstream os;
            os << "'" << fmt_choice << "' is not a valid output format";
            string msg = CNcbiOstrstreamToString(os);
            NCBI_THROW(CInputException, eInvalidInput, msg);
        }
        if (val < 0 || val >= static_cast<int>(eEndValue))
        {
            string msg("Formatting choice is out of range");
            throw std::out_of_range(msg);
        }
        if (m_IsIgBlast && (val != 3 && val != 4 && val != 7))
        {
            string msg("Formatting choice is not valid");
            throw std::out_of_range(msg);
        }
        fmt_type = static_cast<EOutputFormat>(val);
        if ( !(fmt_type == eTabular ||
                fmt_type == eTabularWithComments ||
                fmt_type == eCommaSeparatedValues) )
        {
            custom_fmt_spec.clear();
        }
    }
}


void
CFormattingArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opt)
{
    ParseFormattingString(args, m_OutputFormat, m_CustomOutputFormatSpec);
    m_ShowGis = static_cast<bool>(args[kArgShowGIs]);
    if(m_IsIgBlast)
    {
        m_Html = false;
    }
    else
    {
        m_Html = static_cast<bool>(args[kArgProduceHtml]);
    }
    // Default hitlist size 500, value can be changed if import search strategy is used
    int hitlist_size = opt.GetHitlistSize();

    // To preserve hitlist size in import search strategy > 500,
    // we need to increase the  num_ descriptions and num_alignemtns
    if(hitlist_size > BLAST_HITLIST_SIZE )
    {
        if(!args[kArgNumDescriptions] && !args[kArgNumAlignments] &&
                (m_OutputFormat <= eFlatQueryAnchoredNoIdentities))
        {
            m_NumDescriptions = hitlist_size;
            m_NumAlignments = hitlist_size/ 2;
            return;
        }
    }

    if(m_OutputFormat <= eFlatQueryAnchoredNoIdentities)
    {
        if (args[kArgMaxTargetSequences])
        {
            ERR_POST(Warning << "The parameter -max_target_seqs is ignored for "
                     "output formats, 0,1,2,3. Use -num_descriptions "
                     "and -num_alignments to control output");
        }

        m_NumDescriptions = m_DfltNumDescriptions;
        m_NumAlignments = m_DfltNumAlignments;

        if (args[kArgNumDescriptions])
        {
            m_NumDescriptions = args[kArgNumDescriptions].AsInteger();
        }

        if (args[kArgNumAlignments])
        {
            m_NumAlignments = args[kArgNumAlignments].AsInteger();
        }

        // The If clause is for handling import_search_strategy hitlist size < 500
        // We want to preserve the hitlist size in iss if no formatting input is entered in cmdline
        // If formmating option(s) is entered than the iss hitlist size is overridden.
        if (args[kArgNumDescriptions] || args[kArgNumAlignments])
        {
            hitlist_size = max(m_NumDescriptions, m_NumAlignments);
        }
    }
    else
    {
        if (args[kArgNumDescriptions])
        {
            ERR_POST(Warning << "The parameter -num_descriptions is ignored for "
                     "output formats > 4 . Use -max_target_seqs "
                     "to control output");
        }

        if (args[kArgMaxTargetSequences])
        {
            hitlist_size = args[kArgMaxTargetSequences].AsInteger();
        }
        else if (args[kArgNumAlignments])
        {
            hitlist_size = args[kArgNumAlignments].AsInteger();
        }

        m_NumDescriptions = hitlist_size;
        m_NumAlignments = hitlist_size;
    }

    opt.SetHitlistSize(hitlist_size);

    return;
}

void
CMTArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)


{
    if(m_IsRpsBlast)
    {
        x_SetArgumentDescriptionsRpsBlast(arg_desc);
        return;
    }

    const int kMinValue = static_cast<int>(CThreadable::kMinNumThreads);

    // number of threads
    arg_desc.SetCurrentGroup("Miscellaneous options");
    arg_desc.AddDefaultKey(kArgNumThreads, "int_value",
                           "Number of threads (CPUs) to use in the BLAST search",
                           CArgDescriptions::eInteger,
                           NStr::IntToString(kMinValue));
    arg_desc.SetConstraint(kArgNumThreads,
                           new CArgAllowValuesGreaterThanOrEqual(kMinValue));
    arg_desc.SetDependency(kArgNumThreads,
                           CArgDescriptions::eExcludes,
                           kArgRemote);
    /*
    arg_desc.SetDependency(kArgNumThreads,
                           CArgDescriptions::eExcludes,
                           kArgUseIndex);
    */
    arg_desc.SetCurrentGroup("");
}

void
CMTArgs::x_SetArgumentDescriptionsRpsBlast(CArgDescriptions& arg_desc)
{
    // number of threads
    arg_desc.SetCurrentGroup("Miscellaneous options");
    arg_desc.AddDefaultKey(kArgNumThreads, "int_value",
                           "Number of threads to use in RPS BLAST search:\n "
                           "0 (auto = num of databases)\n "
                           "1 (disable)\n max number of threads = num of databases",
                           CArgDescriptions::eInteger,
                           NStr::IntToString(kDefaultRpsNumThreads));
    arg_desc.SetConstraint(kArgNumThreads,
                           new CArgAllowValuesGreaterThanOrEqual(0));
    arg_desc.SetDependency(kArgNumThreads,
                           CArgDescriptions::eExcludes,
                           kArgRemote);
    /*
    arg_desc.SetDependency(kArgNumThreads,
                           CArgDescriptions::eExcludes,
                           kArgUseIndex);
    */
    arg_desc.SetCurrentGroup("");
}

void
CMTArgs::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& /* opts */)
{
    if(m_IsRpsBlast)
    {
        x_ExtractAlgorithmOptionsRpsBlast(args);
        return;
    }
    if (args.Exist(kArgNumThreads) &&
            args[kArgNumThreads].HasValue())    // could be cancelled by the exclusion in CRemoteArgs
    {
        m_NumThreads = args[kArgNumThreads].AsInteger();

        // This is temporarily ignored (per SB-635)
        if (args.Exist(kArgSubject) && args[kArgSubject].HasValue() &&
                m_NumThreads != CThreadable::kMinNumThreads)
        {
            m_NumThreads = CThreadable::kMinNumThreads;
            LOG_POST(Warning << "'" << kArgNumThreads << "' is currently "
                     << "ignored when '" << kArgSubject << "' is specified.");
        }
    }
}

void
CMTArgs::x_ExtractAlgorithmOptionsRpsBlast(const CArgs& args)
{
    if (args.Exist(kArgNumThreads) &&
            args[kArgNumThreads].HasValue())
    {
        m_NumThreads = args[kArgNumThreads].AsInteger();
    }
}

void
CRemoteArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("Miscellaneous options");
    arg_desc.AddFlag(kArgRemote, "Execute search remotely?", true);

    arg_desc.SetCurrentGroup("");
}

void
CRemoteArgs::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& /* opts */)
{
    if (args.Exist(kArgRemote))
    {
        m_IsRemote = static_cast<bool>(args[kArgRemote]);
    }
}

void
CDebugArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
#if _BLAST_DEBUG
    arg_desc.SetCurrentGroup("Miscellaneous options");
    arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
                     true);
    arg_desc.AddFlag("remote_verbose",
                     "Produce verbose output for remote searches", true);
    arg_desc.AddFlag("use_test_remote_service",
                     "Send remote requests to test servers", true);
    arg_desc.SetCurrentGroup("");
#endif /* _BLAST_DEBUG */
}

void
CDebugArgs::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& /* opts */)
{
#if _BLAST_DEBUG
    m_DebugOutput = static_cast<bool>(args["verbose"]);
    m_RmtDebugOutput = static_cast<bool>(args["remote_verbose"]);
    if (args["use_test_remote_service"])
    {
        IRWRegistry& reg = CNcbiApplication::Instance()->GetConfig();
        reg.Set("BLAST4", DEF_CONN_REG_SECTION "_" REG_CONN_SERVICE_NAME,
                "blast4_test");
    }
#endif /* _BLAST_DEBUG */
}

void
CHspFilteringArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    // culling limit
    arg_desc.SetCurrentGroup("Restrict search or results");
    arg_desc.AddOptionalKey(kArgCullingLimit, "int_value",
                            "If the query range of a hit is enveloped by that of at "
                            "least this many higher-scoring hits, delete the hit",
                            CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgCullingLimit,
                           // best hit algorithm arguments
                           new CArgAllowValuesGreaterThanOrEqual(kDfltArgCullingLimit));

    arg_desc.AddOptionalKey(kArgBestHitOverhang, "float_value",
                            "Best Hit algorithm overhang value "
                            "(recommended value: " +
                            NStr::DoubleToString(kDfltArgBestHitOverhang) +
                            ")",
                            CArgDescriptions::eDouble);
    arg_desc.SetConstraint(kArgBestHitOverhang,
                           new CArgAllowValuesBetween(kBestHit_OverhangMin,
                                   kBestHit_OverhangMax));
    arg_desc.SetDependency(kArgBestHitOverhang,
                           CArgDescriptions::eExcludes,
                           kArgCullingLimit);

    arg_desc.AddOptionalKey(kArgBestHitScoreEdge, "float_value",
                            "Best Hit algorithm score edge value "
                            "(recommended value: " +
                            NStr::DoubleToString(kDfltArgBestHitScoreEdge) +
                            ")",
                            CArgDescriptions::eDouble);
    arg_desc.SetConstraint(kArgBestHitScoreEdge,
                           new CArgAllowValuesBetween(kBestHit_ScoreEdgeMin,
                                   kBestHit_ScoreEdgeMax));
    arg_desc.SetDependency(kArgBestHitScoreEdge,
                           CArgDescriptions::eExcludes,
                           kArgCullingLimit);
    arg_desc.SetCurrentGroup("");
}

void
CHspFilteringArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& opts)
{
    if (args[kArgCullingLimit])
    {
        opts.SetCullingLimit(args[kArgCullingLimit].AsInteger());
    }
    if (args[kArgBestHitOverhang])
    {
        opts.SetBestHitOverhang(args[kArgBestHitOverhang].AsDouble());
    }
    if (args[kArgBestHitScoreEdge])
    {
        opts.SetBestHitScoreEdge(args[kArgBestHitScoreEdge].AsDouble());
    }
}

void
CMbIndexArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("General search options");
    arg_desc.AddOptionalKey(
        kArgUseIndex, "boolean",
        "Use MegaBLAST database index",
        CArgDescriptions::eBoolean );
    arg_desc.AddOptionalKey(
        kArgIndexName, "string",
        "MegaBLAST database index name",
        CArgDescriptions::eString );
    arg_desc.SetCurrentGroup( "" );
}

bool
CMbIndexArgs::HasBeenSet(const CArgs& args)
{
    if ( (args.Exist(kArgUseIndex) && args[kArgUseIndex].HasValue()) ||
            (args.Exist(kArgIndexName) && args[kArgIndexName].HasValue()) )
    {
        return true;
    }
    return false;
}

void
CMbIndexArgs::ExtractAlgorithmOptions(const CArgs& args,
                                      CBlastOptions& opts)
{
    // MB Index does not apply to Blast2Sequences
    if( args.Exist( kArgUseIndex ) &&
            !(args.Exist( kArgSubject ) && args[kArgSubject]))
    {

        bool use_index   = true;
        bool force_index = false;
        bool old_style_index = false;

        if( args[kArgUseIndex] )
        {
            if( args[kArgUseIndex].AsBoolean() ) force_index = true;
            else use_index = false;
        }

        if( args.Exist( kTask ) && args[kTask] &&
                args[kTask].AsString() != "megablast" )
        {
            use_index = false;
        }

        if( use_index )
        {
            string index_name;

            if( args.Exist( kArgIndexName ) && args[kArgIndexName] )
            {
                index_name = args[kArgIndexName].AsString();
                old_style_index = true;
            }
            else if( args.Exist( kArgDb ) && args[kArgDb] )
            {
                index_name = args[kArgDb].AsString();
            }
            else
            {
                NCBI_THROW(CInputException, eInvalidInput,
                           "Can not deduce database index name" );
            }

            opts.SetUseIndex( true, index_name, force_index, old_style_index );
        }
    }
}

void
CStdCmdLineArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("Input query options");

    // query filename
    arg_desc.AddDefaultKey(kArgQuery, "input_file",
                           "Input file name",
                           CArgDescriptions::eInputFile, kDfltArgQuery);

    arg_desc.SetCurrentGroup("General search options");

    // report output file
    arg_desc.AddDefaultKey(kArgOutput, "output_file",
                           "Output file name",
                           CArgDescriptions::eOutputFile, "-");

    arg_desc.SetCurrentGroup("");
}

void
CStdCmdLineArgs::ExtractAlgorithmOptions(const CArgs& args,
        CBlastOptions& /* opt */)
{
    if (args.Exist(kArgQuery) && args[kArgQuery].HasValue() &&
            m_InputStream == NULL)
    {
        m_InputStream = &args[kArgQuery].AsInputFile();
    }
    m_OutputStream = &args[kArgOutput].AsOutputFile();
}

CNcbiIstream&
CStdCmdLineArgs::GetInputStream() const
{
    // programmer must ensure the ExtractAlgorithmOptions method is called
    // before this method is invoked
    if ( !m_InputStream )
    {
        abort();
    }
    return *m_InputStream;
}

CNcbiOstream&
CStdCmdLineArgs::GetOutputStream() const
{
    // programmer must ensure the ExtractAlgorithmOptions method is called
    // before this method is invoked
    _ASSERT(m_OutputStream);
    return *m_OutputStream;
}

void
CStdCmdLineArgs::SetInputStream(CRef<CTmpFile> input_file)
{
    m_QueryTmpInputFile = input_file;
    m_InputStream = &input_file->AsInputFile(CTmpFile::eIfExists_Throw);
}

void
CSearchStrategyArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("Search strategy options");

    arg_desc.AddOptionalKey(kArgInputSearchStrategy,
                            "filename",
                            "Search strategy to use",
                            CArgDescriptions::eInputFile);
    arg_desc.AddOptionalKey(kArgOutputSearchStrategy,
                            "filename",
                            "File name to record the search strategy used",
                            CArgDescriptions::eOutputFile);
    arg_desc.SetDependency(kArgInputSearchStrategy,
                           CArgDescriptions::eExcludes,
                           kArgOutputSearchStrategy);

    arg_desc.SetCurrentGroup("");
}

void
CSearchStrategyArgs::ExtractAlgorithmOptions(const CArgs& /* cmd_line_args */,
        CBlastOptions& /* options */)
{
}

CNcbiIstream*
CSearchStrategyArgs::GetImportStream(const CArgs& args) const
{
    CNcbiIstream* retval = NULL;
    if (args[kArgInputSearchStrategy].HasValue())
    {
        retval = &args[kArgInputSearchStrategy].AsInputFile();
    }
    return retval;
}

CNcbiOstream*
CSearchStrategyArgs::GetExportStream(const CArgs& args) const
{
    CNcbiOstream* retval = NULL;
    if (args[kArgOutputSearchStrategy].HasValue())
    {
        retval = &args[kArgOutputSearchStrategy].AsOutputFile();
    }
    return retval;
}

CBlastAppArgs::CBlastAppArgs()
{
    m_SearchStrategyArgs.Reset(new CSearchStrategyArgs);
    m_Args.push_back(CRef<IBlastCmdLineArgs>(&*m_SearchStrategyArgs));
    m_IsUngapped = false;
}

CArgDescriptions*
CBlastAppArgs::SetCommandLine()
{
    return SetUpCommandLineArguments(m_Args);
}

CRef<CBlastOptionsHandle>
CBlastAppArgs::SetOptions(const CArgs& args)
{
    // We're recovering from a saved strategy or combining
    // CBlastOptions/CBlastOptionsHandle with command line options (in GBench,
    // see GB-1116), so we need to still extract
    // certain options from the command line, include overriding query
    // and/or database
    if (m_OptsHandle.NotEmpty())
    {
        CBlastOptions& opts = m_OptsHandle->SetOptions();
        //opts.DebugDumpText(cerr, "OptionsBeforeLoop", 1);
        const bool mbidxargs_set = CMbIndexArgs::HasBeenSet(args);
        const bool dbargs_set = CBlastDatabaseArgs::HasBeenSet(args);
        NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args)
        {
            if (dynamic_cast<CMbIndexArgs*>(&**arg))
            {
                if (mbidxargs_set)
                    (*arg)->ExtractAlgorithmOptions(args, opts);
            }
            else if (dynamic_cast<CBlastDatabaseArgs*>(&**arg))
            {
                if (dbargs_set)
                    m_BlastDbArgs->ExtractAlgorithmOptions(args, opts);
            }
            else
            {
                (*arg)->ExtractAlgorithmOptions(args, opts);
            }
        }
        m_IsUngapped = !opts.GetGappedMode();
        try
        {
            m_OptsHandle->Validate();
        }
        catch (const CBlastException& e)
        {
            NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
        }
        //opts.DebugDumpText(cerr, "OptionsAfterLoop", 1);
        return m_OptsHandle;
    }

    CBlastOptions::EAPILocality locality =
        (args.Exist(kArgRemote) && args[kArgRemote])
        ? CBlastOptions::eRemote
        : CBlastOptions::eLocal;

    // This is needed as a CRemoteBlast object and its options are instantiated
    // to create the search strategy
    if (GetExportSearchStrategyStream(args) ||
            m_FormattingArgs->ArchiveFormatRequested(args))
    {
        locality = CBlastOptions::eBoth;
    }

    CRef<CBlastOptionsHandle> retval(x_CreateOptionsHandle(locality, args));
    CBlastOptions& opts = retval->SetOptions();
    NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args)
    {
        (*arg)->ExtractAlgorithmOptions(args, opts);
    }

    m_IsUngapped = !opts.GetGappedMode();
    try
    {
        retval->Validate();
    }
    catch (const CBlastException& e)
    {
        NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
    }
    return retval;
}

void CBlastAppArgs::SetTask(const string& task)
{
#if _BLAST_DEBUG
    ThrowIfInvalidTask(task);
#endif
    m_Task.assign(task);
}
/* *********** START ************* */
void
CGpuArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{

    arg_desc.SetCurrentGroup("GPU Options");

	    arg_desc.AddDefaultKey(kArgUseGpu, "int_value",
                           "Use GPU for BLASTX/BLASTP or not, \n"
			   "  0 = not use GPU,\n"
			   "  1...8 = use GPU for BLASTX/BLASTP with the \n"
			   "    GPU card id. \n",
                           CArgDescriptions::eInteger,
                           kDfltArgUseGpu);

    arg_desc.SetConstraint(kArgUseGpu, new CArgAllow_Integers(0, 8) );

    arg_desc.AddDefaultKey(kArgNumGpuThreads, "int_value",
                           "Number of GPU threads per block",
                           CArgDescriptions::eInteger,
                           kDfltArgNumGpuThreads);

    arg_desc.SetConstraint(kArgNumGpuThreads, new CArgAllow_Integers(1,65536) );

    arg_desc.AddDefaultKey(kArgNumGpuBlocks, "int_value",
                           "Number of GPU blocks per grid",
                           CArgDescriptions::eInteger,
                           kDfltArgNumGpuBlocks);

    arg_desc.SetConstraint(kArgNumGpuBlocks, new CArgAllow_Integers(1,65536) );

    arg_desc.AddDefaultKey(kArgGpuMethod, "int_value",
                           "Method to be used\n"
			   "  0 = for GPU-based sequence alignment with all gpu cards in the machine,\n"
                           "  1 = for GPU-based sequence alignment with 1 gpu card (default),\n"
                           "  2...8 = for GPU-based sequence alignment with no. of gpu cards,\n"
                           "  -1 = for GPU database creation, \n",
                           CArgDescriptions::eInteger, kDfltArgGpuMethod);

    arg_desc.SetConstraint(kArgGpuMethod, new CArgAllow_Integers(-1, 8) );

}

void
CGpuArgs::ExtractAlgorithmOptions(const CArgs& args,
                                  CBlastOptions& opts)
{
    if( args[kArgUseGpu] )
    {
        opts.SetUseGpu(args[kArgUseGpu].AsInteger());
    }

    if( args[kArgNumGpuThreads] )
    {
        opts.SetNumGpuThreads(args[kArgNumGpuThreads].AsInteger());
    }

    if( args[kArgNumGpuBlocks] )
    {
        opts.SetNumGpuBlocks(args[kArgNumGpuBlocks].AsInteger());
    }

    //  if( args[kArgNumGpuSequences] ){
    //    opts.SetNumGpuSequences(args[kArgNumGpuSequences].AsInteger());
    //  }

    if( args[kArgGpuMethod] )
    {
        opts.SetGpuMethod(args[kArgGpuMethod].AsInteger());
    }

    if (args.Exist(kArgNumThreads) &&
            args[kArgNumThreads].HasValue())    // could be cancelled by the exclusion in CRemoteArgs
    {
        opts.SetNumCpuThreads(args[kArgNumThreads].AsInteger());
    }

}


/* ********** FINISH ************* */


CArgDescriptions*
SetUpCommandLineArguments(TBlastCmdLineArgs& args)
{
    auto_ptr<CArgDescriptions> retval(new CArgDescriptions);

    // Create the groups so that the ordering is established
    retval->SetCurrentGroup("Input query options");
    retval->SetCurrentGroup("General search options");
    retval->SetCurrentGroup("BLAST database options");
    retval->SetCurrentGroup("BLAST-2-Sequences options");
    retval->SetCurrentGroup("Formatting options");
    retval->SetCurrentGroup("Query filtering options");
    retval->SetCurrentGroup("Restrict search or results");
    retval->SetCurrentGroup("Discontiguous MegaBLAST options");
    retval->SetCurrentGroup("Statistical options");
    retval->SetCurrentGroup("Search strategy options");
    retval->SetCurrentGroup("Extension options");
    /* *********** START ************* */
    retval->SetCurrentGroup("GPU options");
    /* ********** FINISH ************* */
    retval->SetCurrentGroup("");


    NON_CONST_ITERATE(TBlastCmdLineArgs, arg, args)
    {
        (*arg)->SetArgumentDescriptions(*retval);
    }
    return retval.release();
}

void
CBlastAppArgs::x_IssueWarningsForIgnoredOptions(const CArgs& args)
{
    set<string> can_override;
    can_override.insert(kArgQuery);
    can_override.insert(kArgQueryLocation);
    can_override.insert(kArgSubject);
    can_override.insert(kArgSubjectLocation);
    can_override.insert(kArgUseLCaseMasking);
    can_override.insert(kArgDb);
    can_override.insert(kArgDbSize);
    can_override.insert(kArgEntrezQuery);
    can_override.insert(kArgDbSoftMask);
    can_override.insert(kArgDbHardMask);
    can_override.insert(kArgUseIndex);
    can_override.insert(kArgIndexName);
    can_override.insert(kArgStrand);
    can_override.insert(kArgParseDeflines);
    can_override.insert(kArgOutput);
    can_override.insert(kArgOutputFormat);
    can_override.insert(kArgNumDescriptions);
    can_override.insert(kArgNumAlignments);
    can_override.insert(kArgMaxTargetSequences);
    can_override.insert(kArgRemote);
    can_override.insert(kArgNumThreads);
    can_override.insert(kArgInputSearchStrategy);
    can_override.insert(kArgRemote);
    can_override.insert("remote_verbose");
    can_override.insert("verbose");

    // this stores the arguments (and their defaults) that cannot be overriden
    map<string, string> has_defaults;
    EBlastProgramType prog = m_OptsHandle->GetOptions().GetProgramType();
    has_defaults[kArgCompBasedStats] =
        Blast_ProgramIsRpsBlast(prog) ? kDfltArgCompBasedStatsDelta:kDfltArgCompBasedStats;
    // FIX the line below for igblast, and add igblast options
    has_defaults[kArgEvalue] = NStr::DoubleToString(BLAST_EXPECT_VALUE);
    has_defaults[kTask] = m_Task;
    has_defaults[kArgOldStyleIndex] = kDfltArgOldStyleIndex;
    has_defaults[kArgMaxHSPsPerSubject] =
        NStr::IntToString(kDfltArgMaxHSPsPerSubject);
    if (Blast_QueryIsProtein(prog))
    {
        if (NStr::Find(m_Task, "blastp") != NPOS ||
                NStr::Find(m_Task, "psiblast") != NPOS)
        {
            has_defaults[kArgSegFiltering] = kDfltArgNoFiltering;
        }
        else
        {
            has_defaults[kArgSegFiltering] = kDfltArgSegFiltering;
        }
        has_defaults[kArgLookupTableMaskingOnly] =
            kDfltArgLookupTableMaskingOnlyProt;
        has_defaults[kArgGapTrigger] =
            NStr::IntToString((int)BLAST_GAP_TRIGGER_PROT);
    }
    else
    {
        has_defaults[kArgDustFiltering] = kDfltArgDustFiltering;
        has_defaults[kArgLookupTableMaskingOnly] =
            kDfltArgLookupTableMaskingOnlyNucl;
        has_defaults[kArgGapTrigger] =
            NStr::IntToString((int)BLAST_GAP_TRIGGER_NUCL);
    }
    has_defaults[kArgOffDiagonalRange] =
        NStr::IntToString(kDfltOffDiagonalRange);
    has_defaults[kArgMaskLevel] = kDfltArgMaskLevel;
    has_defaults[kArgMaxIntronLength] =
        NStr::IntToString(kDfltArgMaxIntronLength);
    has_defaults[kArgQueryGeneticCode] = NStr::IntToString((int)BLAST_GENETIC_CODE);
    has_defaults[kArgDbGeneticCode] = NStr::IntToString((int)BLAST_GENETIC_CODE);
    // pssm engine/psiblast default options
    has_defaults[kArgPSIPseudocount] =
        NStr::IntToString(PSI_PSEUDO_COUNT_CONST);
    has_defaults[kArgPSIInclusionEThreshold] =
        NStr::DoubleToString(PSI_INCLUSION_ETHRESH);
    has_defaults[kArgPSINumIterations] = NStr::IntToString(1);

    // get arguments, remove the supported ones and warn about those that
    // cannot be overridden.
    typedef vector< CRef<CArgValue> > TArgs;
    TArgs arguments = args.GetAll();
    ITERATE(TArgs, a, arguments)
    {
        const string& arg_name = (*a)->GetName();
        const string& arg_value = (*a)->AsString();
        // if it has a default value, ignore it if it's not different from the
        // default, otherwise, issue a warning
        if (has_defaults.find(arg_name) != has_defaults.end())
        {
            if (has_defaults[arg_name] == arg_value)
            {
                continue;
            }
            else
            {
                if (arg_name == kTask && arg_value == "megablast")
                {
                    // No need to issue warning here, as it's OK to change this
                    continue;
                }
                LOG_POST(Warning << arg_name << " cannot be overridden when "
                         "using a search strategy");
            }
        }
        // if the argument cannot be overridden, issue a warning
        if (can_override.find(arg_name) == can_override.end())
        {
            LOG_POST(Warning << arg_name << " cannot be overridden when "
                     "using a search strategy");
        }
    }
}

CRef<CBlastOptionsHandle>
CBlastAppArgs::SetOptionsForSavedStrategy(const CArgs& args)
{
    if(m_OptsHandle.Empty())
    {
        NCBI_THROW(CInputException, eInvalidInput, "Empty Blast Options Handle");
    }

    // We're recovering from a saved strategy, so we need to still extract
    // certain options from the command line, include overriding query
    // and/or database
    CBlastOptions& opts = m_OptsHandle->SetOptions();
    // invoke ExtractAlgorithmOptions on certain argument classes, i.e.: those
    // that should have their arguments overriden
    m_QueryOptsArgs->ExtractAlgorithmOptions(args, opts);
    m_StdCmdLineArgs->ExtractAlgorithmOptions(args, opts);
    m_RemoteArgs->ExtractAlgorithmOptions(args, opts);
    m_DebugArgs->ExtractAlgorithmOptions(args, opts);
    m_FormattingArgs->ExtractAlgorithmOptions(args, opts);
    m_MTArgs->ExtractAlgorithmOptions(args, opts);
    if (CBlastDatabaseArgs::HasBeenSet(args))
    {
        m_BlastDbArgs->ExtractAlgorithmOptions(args, opts);
    }
    if (CMbIndexArgs::HasBeenSet(args))
    {
        NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args)
        {
            if (dynamic_cast<CMbIndexArgs*>(arg->GetPointer()) != NULL)
            {
                (*arg)->ExtractAlgorithmOptions(args, opts);
            }
        }
    }
    m_IsUngapped = !opts.GetGappedMode();
    x_IssueWarningsForIgnoredOptions(args);
    try
    {
        m_OptsHandle->Validate();
    }
    catch (const CBlastException& e)
    {
        NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
    }
    return m_OptsHandle;
}

END_SCOPE(blast)
END_NCBI_SCOPE
