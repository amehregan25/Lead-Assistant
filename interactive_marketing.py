"""
interactive_marketing.py
==========================

This script provides an interactive interface for generating tailored
recommendation reports based on Nuclera’s case study library.  Users can
input a lead’s first and last name, job title, company, workflow (soluble,
membrane or antibody), protein targets and therapeutic areas, plus an
optional LinkedIn summary.  The script infers an appropriate market
segment (pharma, biotech SME or academia) and persona (decision maker,
scientist or core facility manager) based on the company, title and
summary.  It maps the workflow to one of several scientific categories
(membrane, antibody, kinase, transcription factor or generic) to ensure
that recommendations align with the lead’s work.

Rather than composing emails, this module focuses on sharing relevant data
and recommending case studies.  For a given input, the generated report
returns:

* Whether the lead’s company already has an installed eProtein
  Discovery™ system (current customer) and the location of that unit.
* A list of case studies that directly match the lead’s protein class or
  therapeutic area, along with concise summaries.  If no match is found
  the report falls back to generic workflow recommendations.
* Suggested collateral items (flyers, application notes, slides) for the
  inferred workflow.
* A record of any marketing campaigns the lead has previously been a part
  of, based on a local CSV file mapping first and last names to
  campaign names and first‑association dates.

This module can be run in two ways:

1. As a standalone script: execute ``python interactive_marketing.py`` and
   follow the prompts to enter lead information.  The program prints the
   recommendation report to the console.
2. As an importable library: call ``generate_recommendation_report(
   title, company, topic_input, protein_targets, therapeutic_areas,
   linkedin_summary, first_name, last_name)`` from your own code to
   receive a dictionary containing the mapped topic, segment, persona,
   current customer status and location, recommended case studies,
   suggested resources and prior campaign history.  In this context
   ``topic_input`` refers to the workflow (soluble, membrane or antibody)
   and can be left blank to default to a generic workflow.

Earlier versions of this file included an email campaign generator.  That
functionality has been removed in favour of data sharing and case study
recommendation.  If email generation is required in future, a new
function can be added without affecting the existing recommendation logic.

"""

from __future__ import annotations

from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------
# Normalisation and configuration
#
# To classify companies as “big pharma” we maintain a small set of
# well‑known pharmaceutical corporations.  Leads from these companies
# are assigned to the pharma segment even if their titles don’t
# explicitly mention “pharma”.

def _norm(text: str) -> str:
    """Lowercase and strip whitespace from text for simple matching."""
    return (text or "").lower().strip()


BIG_PHARMA_COMPANIES = {
    "pfizer",
    "novartis",
    "roche",
    "merck",
    "msd",
    "bristol-myers squibb",
    "bristol myers",
    "bms",
    "johnson & johnson",
    "janssen",
    "gsk",
    "glaxosmithkline",
    "sanofi",
    "astrazeneca",
    "eli lilly",
    "abbvie",
    "amgen",
    "bayer",
    "takeda",
    "boehringer ingelheim",
    "novo nordisk",
    "regeneron",
    "biogen",
    "vertex",
}

# ---------------------------------------------------------------------------
# Installed base
#
# For current accounts, the script can tailor emails to encourage leads to
# leverage an existing eProtein Discovery™ unit at their organisation.  The
# following dictionary maps account names (lower‑cased) to the location of
# the installed instrument.  If a lead’s company name contains one of these
# account names (case‑insensitive), the email copy will invite them to use
# the unit already available on site.
INSTALL_BASE = {
    "astrazeneca": "Cambridge",
    "diamond": "Oxford",
    "ribbon": "Germany",
    "embl": "Heidelberg",
    "institut pasteur": "Paris",
    "vib": "Ghent",
    "ku leuven": "Leuven",
    "harvard": "Longwood",
    "broad institute": "Harvard Medical School",
    "wyss institute": "Harvard Medical School",
    "mit": "Harvard Medical School",
    "bms": "Princeton",
    "northwestern": "Chicago",
    "stowers": "Kansas City",
    "children's mercy hospital": "Kansas City",
    "monterosa": "Boston",
    "pfizer": "Groton",
    "ucsf": "San Francisco",
    "eli lilly": "Colorado",
}

# ---------------------------------------------------------------------------
# Campaign history loading
#
# The attached CSV file contains campaign membership records with columns:
# Campaign ID, Campaign Name, Member Type, Member First Associated Date,
# First Name, Last Name and Related Record ID.  We load this file once
# into a dictionary keyed by (first_name.lower(), last_name.lower()) to
# quickly look up any campaigns a lead has participated in.  The file is
# expected to be located in the shared folder under the name
# ``report1765405239456.csv``.

import csv
from pathlib import Path

BASE_DIR = Path("campaign_history.csv").resolve().parent

_campaign_data_loaded = False
_campaign_data: Dict[Tuple[str, str], List[Tuple[str, str]]] = {}

def _load_campaign_data() -> None:
    """
    Read the campaign history CSV and populate the global dictionary.  Each
    entry is keyed by (first_name_lower, last_name_lower) and maps to a
    list of (campaign_name, first_associated_date) tuples.
    """
    global _campaign_data_loaded, _campaign_data
    if _campaign_data_loaded:
        return
    csv_path = BASE_DIR / "campaign_history.csv"
    if not csv_path.exists():
        # No campaign data available
        _campaign_data_loaded = True
        return
    try:
        with csv_path.open(newline="", encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f)
            for row in reader:
                fn = row.get("First Name", "").strip().lower()
                ln = row.get("Last Name", "").strip().lower()
                if not fn or not ln:
                    continue
                campaign_name = row.get("Campaign Name", "").strip()
                date = row.get("Member First Associated Date", "").strip()
                key = (fn, ln)
                _campaign_data.setdefault(key, []).append((campaign_name, date))
        _campaign_data_loaded = True
    except Exception:
        # If reading fails for any reason, mark as loaded to avoid repeated attempts
        _campaign_data_loaded = True

def get_campaign_history(first_name: str, last_name: str) -> List[Dict[str, str]]:
    """
    Return a list of campaigns that match the given first and last name.

    Parameters
    ----------
    first_name : str
        The lead’s first name.
    last_name : str
        The lead’s last name.

    Returns
    -------
    List[Dict[str, str]]
        Each entry contains ``campaign_name`` and ``date``.  If no campaigns
        are found, returns an empty list.
    """
    _load_campaign_data()
    fn = first_name.strip().lower() if first_name else ""
    ln = last_name.strip().lower() if last_name else ""
    if not fn or not ln:
        return []
    matches = _campaign_data.get((fn, ln), [])
    return [{"campaign_name": name, "date": date} for name, date in matches]

def get_install_location(company: str) -> str | None:
    """
    Check whether the given company corresponds to an installed eProtein
    Discovery™ instrument.

    Parameters
    ----------
    company : str
        Name of the lead’s company or organisation.

    Returns
    -------
    Optional[str]
        The location of the installed system if the company name matches an
        account in INSTALL_BASE (case‑insensitive, partial match).  Returns
        None if no match is found.
    """
    if not company:
        return None
    company_n = _norm(company)
    for account, location in INSTALL_BASE.items():
        account_n = _norm(account)
        # Consider both directions for partial match: the lead’s company may
        # include the account name or vice versa (e.g. “AstraZeneca plc”).
        if account_n in company_n or company_n in account_n:
            return location
    return None


# ---------------------------------------------------------------------------
# Topic mapping
#
# The following lists define keywords that map a free‑form “topic”
# string to one of the internal categories used by this script.  A lead
# studying GPCRs or ion channels should trigger the “membrane” topic,
# while a lead working on antibody discovery triggers the “antibody”
# topic.  Kinase, transcription factor and generic categories are
# similarly inferred.

MEMBRANE_KEYWORDS = [
    "membrane", "gpcr", "g protein", "ion channel", "channel", "transporter",
    "abc transporter", "nanodisc", "cryo", "cryo-em", "proteoliposome",
]

ANTIBODY_KEYWORDS = [
    "antibody", "nanobody", "vhh", "scfv", "single-chain", "mini-binder",
    "minibinder", "b cell", "b-cell", "display", "phage display",
    "yeast display", "spycatcher", "spytag", "spybli", "binder",
]

KINASE_KEYWORDS = [
    "kinase", "btk", "tyrosine kinase", "serine/threonine kinase",
    "kinome", "phosphorylation", "signal transduction",
]

TF_KEYWORDS = [
    "transcription factor", "tf", "chromatin", "gata", "jun", "fos",
    "nrf2", "nfe2l2", "dna binding", "promoter", "enhancer",
]


def map_topic(topic_input: str) -> str:
    """
    Map a free‑form topic description to an internal category.

    Parameters
    ----------
    topic_input : str
        The user‑provided description of the research topic.  This can
        include terms like “GPCR”, “antibody”, “kinase”, etc.

    Returns
    -------
    str
        One of ``"membrane"``, ``"antibody"``, ``"kinase"``, ``"tf"`` or
        ``"generic"``.
    """
    text = _norm(topic_input)
    # Check in priority order: antibody → membrane → kinase → tf
    # Treat both antibody-specific keywords and broad biologic terms as antibody topics.
    for kw in ANTIBODY_KEYWORDS:
        if kw in text:
            return "antibody"
    # Recognise generic biologic/biologics terms as antibody workflows too
    if "biologic" in text or "biologics" in text:
        return "antibody"
    for kw in MEMBRANE_KEYWORDS:
        if kw in text:
            return "membrane"
    for kw in KINASE_KEYWORDS:
        if kw in text:
            return "kinase"
    for kw in TF_KEYWORDS:
        if kw in text:
            return "tf"
    return "generic"


# ---------------------------------------------------------------------------
# Segment and persona inference
#
# These functions use the lead’s company name and title to infer an
# appropriate market segment (pharma, biotech_sme, academia) and
# persona.  Decision makers are identified by keywords such as
# “director” or “vice president”, while core facility managers are
# recognised by phrases like “core facility”.

DECISION_MAKER_KEYWORDS = [
    "vp", "vice president", "director", "head of", "global head",
    "site head", "chief", "cso", "cto", "ceo", "founder", "co-founder",
    "president", "executive", "senior director", "senior manager",
    # scientist titles will be handled by SCIENTIST_KEYWORDS below
]

# Scientist persona keywords.  If a lead’s title contains one of these terms and the
# company is not a recognised big pharma and the lead is not a core facility manager,
# we classify them as a scientist in the biotech SME segment.  Titles containing these
# keywords should be treated as scientists rather than decision makers.
SCIENTIST_KEYWORDS = [
    "lead scientist", "principal scientist", "scientist", "r&d specialist",
    "associate scientist", "research specialist", "researcher",
]

CORE_FACILITY_KEYWORDS = [
    "core facility", "core lab", "shared resource", "platform lead",
    "facility director", "facility manager",
]

ACADEMIA_HINTS = [
    "university", "institute", "college", "hospital", "school",
    "research center", "core facility", "core lab",
]


def infer_segment_persona(company: str, title: str, summary: str = "") -> Tuple[str, str]:
    """
    Infer the market segment and persona from company and title.

    Parameters
    ----------
    company : str
        Name of the lead’s company or organisation.
    title : str
        The lead’s job title.

    Returns
    -------
    Tuple[str, str]
        A tuple ``(segment, persona)`` where ``segment`` is one of
        ``"pharma"``, ``"biotech_sme"`` or ``"academia"`` and ``persona``
        is ``"decision_maker"``, ``"scientist"`` or ``"core_facility"``.

    The logic prioritises recognising scientists working in non‑pharma organisations.
    If the company is not a recognised big pharma and the title contains one of
    ``SCIENTIST_KEYWORDS``, the persona is set to "scientist" and the segment to
    ``"biotech_sme"`` unless academia hints apply.  Decision maker keywords are
    applied only when scientist keywords are not present.
    """
    company_n = _norm(company)
    title_n = _norm(title)
    summary_n = _norm(summary)
    # Combine company, title and LinkedIn summary for keyword searches.  This allows
    # detection of academia hints, decision maker keywords and scientist keywords
    # from the free‑form summary text as well as the title and company name.
    text = f"{company_n} {title_n} {summary_n}"

    # Detect whether the company is a recognised big pharma
    big_pharma_match = False
    for name in BIG_PHARMA_COMPANIES:
        if name in company_n:
            big_pharma_match = True
            break

    # Default segment assignment
    if big_pharma_match:
        segment = "pharma"
    else:
        # Look for academia hints first
        segment = "biotech_sme"
        for hint in ACADEMIA_HINTS:
            if hint in text:
                segment = "academia"
                break

    # Default persona assignment
    persona = "scientist"

    # Core facility overrides everything else
    for kw in CORE_FACILITY_KEYWORDS:
        if kw in text:
            persona = "core_facility"
            break
    else:
        # If the title contains any scientist keywords and the company is not big pharma
        # or academia, treat as scientist and assign segment to biotech_sme
        scientist_match = any(kw in title_n for kw in SCIENTIST_KEYWORDS)
        if scientist_match and not big_pharma_match and segment != "academia":
            persona = "scientist"
            segment = "biotech_sme"
        else:
            # Otherwise, check for decision maker keywords
            for kw in DECISION_MAKER_KEYWORDS:
                if kw in text:
                    persona = "decision_maker"
                    break

            # If not big pharma and not academia, and not already scientist, consider
            # generic pharma strings for segment classification
            if not big_pharma_match and segment == "biotech_sme":
                if "pharma" in text or "pharmaceutical" in text or "biopharma" in text:
                    segment = "pharma"

    return segment, persona


# ---------------------------------------------------------------------------
# Citation texts
#
# Each scientific topic has an associated explanatory paragraph with
# citations.  These paragraphs summarise key insights from Nuclera’s
# public materials and include tether IDs for reproducibility.  When
# composing personalised emails, the script appends one of these
# paragraphs to the body so that recipients receive tangible data points
# along with contextualised benefits.

# Descriptions of the platform's capabilities for each topic (citations removed).  These
# paragraphs summarise key insights from Nuclera’s public materials without
# including tether IDs.  They are appended to the email body to provide
# context and data without explicit citations.
DESCRIPTION_TEXTS: Dict[str, str] = {
    "membrane": (
        "Membrane proteins represent a majority of drug targets but are notoriously difficult to express. "
        "Nuclera’s membrane protein workflow screens multiple DNA constructs and expression conditions within 24 hours and yields high‑yield, functional membrane proteins in 48 hours, validated by ATPase assays and cryo‑EM analysis. "
        "Automated nanodroplet screening across dozens of data points maps optimal nanodisc compositions and stabilising additives, delivering assay‑ready membrane proteins in 48 hours."
    ),
    "antibody": (
        "Obtaining reliable binding kinetics is often limited by purification bottlenecks. "
        "SpyBLI couples SpyCatcher003–SpyTag003 to capture bait proteins directly from crude cell‑free lysates, enabling kinetic analysis without purification. "
        "When combined with Nuclera’s eProtein Discovery, high‑quality binding kinetics can be obtained within 24 hours, compressing the clone–purify–measure cycle. "
        "The platform also supports rapid production of functional proteins for kinases and other targets, providing purified material and kinetic data within five days."
    ),
    "kinase": (
        "Our eProtein Discovery system accelerates kinase research by rapidly producing high‑yield Bruton’s Tyrosine Kinase and other kinases. "
        "By screening multiple variants and conditions using digital microfluidics, the system delivers purified, functional kinase proteins in 48 hours and integrates with Biacore™ SPR to provide kinetic characterisation from DNA to data in five days."
    ),
    "tf": (
        "Transcription factors such as GATA‑1 and c‑Jun are notoriously difficult to express due to intrinsic disorder and aggregation. "
        "Nuclera’s platform uses solubility tag screening and customised cell‑free conditions to express and purify these factors within 48 hours, and DNA‑binding assays confirm their activity."
    ),
    "generic": (
        "Traditional protein workflows rely on weeks of trial‑and‑error across constructs, tags, hosts and buffers, delaying progress. "
        "Nuclera’s eProtein Discovery screens up to 24 constructs across eight conditions—yielding 192 expression and 30 purification data points—in just 24 hours and uses automated, cell‑free multiplex screening to deliver assay‑ready proteins in 48 hours. "
        "By customising cell‑free blends with user‑specified additives such as reducing agents, detergents, cofactors and nanodiscs, researchers can optimise folding and yields for difficult proteins. "
        "For disulfide‑rich or cofactor‑dependent proteins, the platform’s 192‑condition screen and additives like PDI and GSSG guide in‑vivo expression and produce purified proteins within 48 hours."
    ),
}

# Suggested collateral for each topic.  These are the names of case studies,
# application notes or flyers from Nuclera’s resource library that the user
# might attach to the emails.  When generating a campaign, the script will
# list the relevant items for the selected topic so that marketers can
# include them as supporting materials.
RESOURCE_SUGGESTIONS: Dict[str, List[str]] = {
    "membrane": [
        "Functional membrane proteins ready in 48 Hours flyer",
        "MsbA application note",
        "Rapidly express optimised membrane proteins in 48 hours webinar",
    ],
    "antibody": [
        "Accelerating binding kinetics with SpyBLI & eProtein Discovery publication",
        "Accelerated kinase drug discovery (BTK) application note",
        "Rapid cell‑free expression of VEGF case study",
        "Antibody services slides",
    ],
    "kinase": [
        "Accelerated kinase drug discovery (BTK) application note",
        "Accelerating binding kinetics with SpyBLI & eProtein Discovery publication",
    ],
    "tf": [
        "Transcription factor poster (GATA‑1 & c‑Jun)",
        "Guiding in vivo expression application note",
    ],
    "generic": [
        "Protein in 48 hours flyer (general workflow)",
        "Soluble protein workflow flyer",
        "Custom additives manual",
        "Guiding in vivo expression application note",
        "Phytoglobin case study",
    ],
}

# Detailed mapping of case studies to protein classes, therapeutic areas and
# concise summaries.  Each entry associates a case with a set of
# protein‑class keywords, a set of therapeutic‑area keywords and a summary
# describing the result of the experiment.  When generating
# recommendations, the script compares the lead’s protein targets and
# therapeutic areas against these keywords to identify relevant case
# studies.
CASE_STUDIES_INFO: Dict[str, Dict[str, object]] = {
    # Soluble case studies
    "vegf": {
        "protein_class": ["vegf", "growth factor", "vascular endothelial growth factor"],
        "therapeutic_areas": ["cancer", "oncology", "eye", "retinopathy", "macular", "inflammation", "inflammatory"],
        "summary": (
        "The VEGF case study demonstrates that adding redox additives such as TRXB1 and PDI/GSSG and screening solubility tags enables correct disulfide formation.  "
        "Predicted small‑scale yields (~26 µg per 200 µL reaction) matched scaled‑up production (~627 µg per mL), delivering about 52 µg of >95 % pure VEGF with EC₅₀ ≈ 12.5 ng/mL in a dimerisation assay.  "
        "These results show how Nuclera’s screen accurately forecasts scale‑up and yields active VEGF in 48 hours."
        ),
    },
    "hrp": {
        "protein_class": ["hrp", "horseradish peroxidase", "peroxidase", "oxidoreductase"],
        "therapeutic_areas": ["diagnostic", "diagnostics", "immunoassay", "cancer"],
        "summary": (
        "In the horseradish peroxidase (HRP) case, co‑expression with hemin, calcium ions, oxidising additives (PDI/GSSG) and a TRX solubility tag increased expression more than 30‑fold compared with tag‑free conditions.  "
        "Removing the redox additives or hemin decreased activity by roughly 70 %, emphasising the need for proper cofactors.  "
        "The screen correctly predicted improved yields, and scale‑up provided milligram quantities of active HRP within 48 hours for diagnostic and immunoassay applications."
        ),
    },
    "mmp1": {
        "protein_class": ["mmp", "mmp1", "matrix metalloproteinase", "metalloprotease"],
        "therapeutic_areas": ["cancer", "oncology", "metastasis", "collagen", "fibrosis"],
        "summary": (
        "The MMP‑1 (matrix metalloproteinase 1) case compared full‑length and truncated constructs while supplementing zinc and calcium ions together with redox additives such as PDI/GSSG.  "
        "Deleting the hemopexin repeats increased predicted small‑scale yield from 16 µg to about 66 µg per 200 µL reaction and improved stability.  "
        "Scale‑up experiments confirmed higher expression levels for the truncated variant and enabled production of active, soluble MMP‑1 within 48 hours for studies on collagen degradation and cancer metastasis."
        ),
    },
    "cd19": {
        "protein_class": ["cd19", "b cell", "b-cell", "b lymphocyte antigen"],
        "therapeutic_areas": ["b cell malignancy", "leukemia", "lymphoma", "car-t", "immunotherapy"],
        "summary": (
        "In the CD19 case, full‑length and truncated variants of the B‑lymphocyte antigen were co‑expressed with chaperones including PDI/GSSG and 3C protease.  "
        "The optimised conditions boosted expression of the extracellular domain (ECD) to roughly 183 µg mL⁻¹ and the intracellular domain (ICD) to about 65 µg mL⁻¹ during scale‑up.  "
        "Purified ECD at high purity was obtained within 48 hours and used to guide subsequent in‑vivo expression, supporting CAR‑T and antibody development programmes."
        ),
    },
    "tgfb1": {
        "protein_class": ["tgf", "tgf-β", "tgf-beta", "transforming growth factor", "tgfb1"],
        "therapeutic_areas": ["cancer", "oncology", "connective tissue", "skeletal", "bone", "musculoskeletal"],
        "summary": (
        "The TGF‑β1 case explored multiple constructs and solubility tags under oxidising conditions to correctly form the cytokine’s three disulfide bonds.  "
        "The variant comprising residues 279–390 produced the highest expression when co‑expressed with oxidising additives such as PDI/GSSG.  "
        "Scale‑up experiments using CUSF and FH8 solubility tags delivered 30 µg and 26 µg of purified TGF‑β1, respectively, with predicted yields of 32.53 µg and 31.94 µg matching experimental values.  "
        "Purity exceeded 83 % and 94 % for the two tags, and active TGF‑β1 was available in hand within 48 hours for downstream studies."
        ),
    },
    "mdm2": {
        "protein_class": ["mdm2", "e3 ligase", "e3 ubiquitin ligase", "p53 ligase"],
        "therapeutic_areas": ["cancer", "oncology", "protac", "targeted degradation"],
        "summary": (
        "In the MDM2 case, AlphaFold‑guided construct design and solubility tag screening were used to overcome the protein’s aggregation‑prone regions.  "
        "Truncating non‑structured loops improved expression by roughly 5.5‑fold and, when paired with chaperones such as 3C protease, Ca2+, Zn2+, TRX‑B1 and DnaK, delivered soluble protein.  "
        "Scale‑up reactions yielded about 25 µg mL⁻¹ of 95 % pure MDM2 within 48 hours.  "
        "An auto‑ubiquitination assay confirmed functionality, providing active material for PROTAC development and p53 stabilisation studies."
        ),
    },
    "parkin": {
        "protein_class": ["parkin", "e3 ligase", "e3 ubiquitin ligase", "parkin e3 ligase"],
        "therapeutic_areas": ["neurodegeneration", "parkinson", "parkinson's", "mitochondrial quality control"],
        "summary": (
        "For Parkin, an E3 ligase involved in mitochondrial quality control, AlphaFold models were used to design constructs with loop deletions.  "
        "Removing an entire loop and appending a ZZ solubility tag substantially improved yields.  "
        "The optimised construct produced around 40 µg mL⁻¹ of purified, soluble Parkin within 48 hours, supplying enough active protein for biochemical assays relevant to early‑onset Parkinson’s disease."
        ),
    },
    "vhl": {
        "protein_class": ["vhl", "von hippel‑lindau", "e3 ligase", "e3 ubiquitin ligase"],
        "therapeutic_areas": ["vhl disease", "renal cancer", "oncology", "protac", "targeted degradation"],
        "summary": (
        "VHL is a substrate‑recognition component of an E3 ligase complex that targets hypoxia‑inducible factor (HIF) for degradation but is prone to instability and aggregation.  "
        "By engineering constructs guided by AlphaFold and screening multiple solubility tags and additives, the platform stabilised VHL for expression.  "
        "The resulting optimised variant yielded functional VHL protein within 48 hours, enabling investigation of von Hippel–Lindau disease, renal cancer and the design of PROTACs that hijack the VHL–HIF interaction."
        ),
    },
    # Membrane case studies
    "ffar1": {
        "protein_class": ["ffar1", "gpr40", "gpcr", "free fatty acid receptor 1","gpcr","gipr","a2ar","orphan gpcr"],
        "therapeutic_areas": ["type 2 diabetes", "diabetes", "obesity", "cancer", "alzheimer", "neurodegenerative"],
        "summary": (
        "The FFAR1 (GPR40) case used nanodisc and detergent screening to identify conditions that stabilised this G‑protein‑coupled receptor.  "
        "The optimised construct produced active receptor at roughly 0.12 mg mL⁻¹ with purity exceeding 90 %.  "
        "Negative‑stain electron microscopy revealed homogeneous particles, demonstrating suitability for structural studies.  "
        "These results were achieved in 48 hours, supporting research into type 2 diabetes, obesity, cancer and neurodegenerative diseases."
        ),
    },
    "msba": {
        "protein_class": ["msba", "abc transporter", "lipid flippase", "abc transporter", "lipid a‑core flippase"],
        "therapeutic_areas": ["gram‑negative", "antibiotic", "antibacterial", "drug resistance"],
        "summary": (
        "For MsbA, an essential lipid A‑core flippase in Gram‑negative bacteria, screening multiple DNA constructs and nanodisc compositions improved predicted yields by roughly five‑fold.  "
        "The optimal conditions produced highly pure MsbA in 48 hours; ATPase assays confirmed activity and cryo‑EM analysis verified correct folding and dimeric assembly.  "
        "These results highlight how nanodisc optimisation enables rapid production of functional ABC transporters for antibiotic discovery."
        ),
    },
    "zmpste24": {
        "protein_class": ["zmpste24", "membrane metalloprotease", "zmp ste24"],
        "therapeutic_areas": ["premature aging", "progeria", "cardiovascular", "nuclear envelope", "laminopathy"],
        "summary": (
        "ZMPSTE24 is a seven‑transmembrane metalloprotease implicated in premature aging disorders, cardiovascular disease and defects of the nuclear envelope.  "
        "To overcome the challenges posed by its helical barrel, the screen varied nanodisc compositions and detergents, identifying conditions that stabilised the protein.  "
        "This approach yielded correctly folded, functional ZMPSTE24 within 48 hours for mechanistic studies of laminopathies and progeroid syndromes."
        ),
    },
    "bacteriorhodopsin": {
        "protein_class": ["bacteriorhodopsin", "br", "light‑driven proton pump", "photoreceptor", "rhodopsin"],
        "therapeutic_areas": ["bioelectronics", "optogenetics", "biomimetics", "energy harvesting"],
        "summary": (
        "Bacteriorhodopsin, a seven‑helix photoreceptor that pumps protons upon absorbing light, requires specific lipids and a retinal chromophore for stability.  "
        "The case study screened nanodisc lipids and retinal concentrations to maximise activity, yielding active BR within 48 hours.  "
        "This rapid production enables applications in bioelectronics, biomimetics and optogenetics that typically rely on more complex expression systems."
        ),
    },
    "dhodh": {
        "protein_class": ["dhodh", "dihydroorotate dehydrogenase", "mitochondrial enzyme"],
        "therapeutic_areas": ["cancer", "oncology", "autoimmune", "pyrimidine synthesis", "metabolic"],
        "summary": (
        "Dihydroorotate dehydrogenase (DHODH) is an inner mitochondrial membrane enzyme essential for de novo pyrimidine synthesis and a target for cancer and autoimmune therapies.  "
        "The case study varied cofactor levels (FMN) and detergents such as Brij®‑35 to solubilise and activate this otherwise insoluble enzyme.  "
        "Soluble, functional DHODH was produced in 48 hours, demonstrating that appropriate cofactor supplementation can rescue challenging membrane enzymes and underscoring that roughly one‑third of human proteins require cofactors."
        ),
    },
    # Kinase and antibody case studies
    "btk": {
        "protein_class": ["btk", "bruton", "tyrosine kinase", "kinase"],
        "therapeutic_areas": ["cancer", "oncology", "autoimmune", "inflammatory"],
        "summary": (
        "In the BTK case, multiple DNA constructs and expression conditions were screened to maximise yield and solubility for Bruton’s tyrosine kinase.  "
        "The top‑performing construct yielded purified BTK within 48 hours; subsequent Biacore™ SPR measurements provided kinetic data in just five days from DNA to decision.  "
        "This workflow accelerates kinase projects by producing functional protein and binding data far more quickly than traditional cell‑based methods."
        ),
    },
    # Transcription factor case studies
    "gata": {
        "protein_class": ["gata", "gata‑1", "c‑jun", "jun", "transcription factor", "tf"],
        "therapeutic_areas": ["leukemia", "blood", "cancer", "oncology", "immunology"],
        "summary": (
        "Transcription factors GATA‑1 and c‑Jun are intrinsically disordered and notoriously difficult to express.  "
        "Through solubility tag screening and redox additive optimisation, the platform produced soluble and active GATA‑1 and c‑Jun within 48 hours.  "
        "DNA‑binding assays confirmed their functionality, demonstrating that the workflow can rapidly generate intrinsically disordered TFs for structural and functional studies."
        ),
    },
}

# Human‑readable titles for each case study.  These descriptive names are used
# in recommendation reports so that users can easily recognise the subject of
# each case.  The keys correspond to entries in ``CASE_STUDIES_INFO``.
CASE_TITLES: Dict[str, str] = {
    "vegf": "VEGF (Vascular Endothelial Growth Factor) case study",
    "hrp": "Horseradish Peroxidase (HRP) case study",
    "mmp1": "MMP‑1 (Matrix Metalloproteinase 1) case study",
    "cd19": "CD19 (B‑lymphocyte antigen) case study",
    "tgfb1": "TGF‑β1 (Transforming Growth Factor beta 1) case study",
    "mdm2": "MDM2 E3 ubiquitin ligase case study",
    "parkin": "Parkin E3 ubiquitin ligase case study",
    "vhl": "VHL (Von Hippel‑Lindau) case study",
    "ffar1": "FFAR1/GPR40 GPCR case study",
    "msba": "MsbA ABC transporter case study",
    "zmpste24": "ZMPSTE24 membrane metalloprotease case study",
    "bacteriorhodopsin": "Bacteriorhodopsin photoreceptor case study",
    "dhodh": "DHODH (Dihydroorotate Dehydrogenase) case study",
    "btk": "BTK (Bruton’s Tyrosine Kinase) case study",
    "gata": "Transcription factor case study (GATA‑1 & c‑Jun)",
}

# Additional description for Nuclera antibody services.  This text summarises our
# recombinant antibody screening service, which uses E. coli cell‑free protein
# synthesis and digital microfluidics to express and screen 96 antibody variants
# in parallel.  The workflow delivers binding data in about 24 hours, yields
# roughly seven times more variants at the same cost as traditional methods
# and produces decision‑grade IgG1 outputs ready for CHO scale‑up.
ANTIBODY_SERVICES_DESCRIPTION = (
    "Traditional secondary screening for antibodies is slow and costly because each clone must be expressed, purified and analysed in series.  "
    "Nuclera’s recombinant antibody services integrate engineered E. coli cell‑free expression with digital microfluidics to express and purify up to 96 antibody variants overnight.  "
    "On‑cartridge fluorescence assays provide binding data within 24 hours, enabling you to rank clones quickly, and the yields are sufficient to produce IgG1 formats ready for CHO scale‑up.  "
    "Because the service screens many variants simultaneously, you can explore seven times more variants at the same cost, reducing attrition early and accelerating the path to lead molecules."
)

def infer_research_statement(protein_targets: str) -> str:
    """
    Generate a brief research statement based on the first protein listed.

    The function recognises a few common classes of proteins and returns
    comments on expression or biological context.  If the protein cannot
    be classified, it returns a generic question to prompt discussion.

    Parameters
    ----------
    protein_targets : str
        Comma‑separated string of protein targets.

    Returns
    -------
    str
        A sentence summarising relevant research or posing a question.
    """
    if not protein_targets:
        return ""
    first = protein_targets.split(",")[0].strip()
    first_lower = first.lower()
    # GPCRs and related membrane receptors
    gpcr_keywords = ["gpcr", "gipr", "gpr", "ffar1", "a2a", "a2ar", "glp1r", "glp-1r", "gper", "cb1", "ccr5", "cxcr4"]
    # Ion channels
    ion_channel_keywords = ["channel", "nav", "kv", "kir", "trp", "scn", "cacna", "kcn", "hc", "p2x", "herg", "ampa", "nmda"]
    # Kinases
    kinase_keywords = ["kinase", "btk", "jak", "braf", "flt3", "src", "mapk", "jak2", "pi3k", "akt"]
    # Transcription factors
    tf_keywords = ["transcription factor", "gata", "jun", "c-jun", "fos", "nrf", "nfe2l2", "tf"]
    # Disulfide‑rich / cytokines
    disulfide_keywords = ["vegf", "epo", "interferon", "il-2", "il2", "il-4", "il-10", "tpa"]
    # Peroxidases / oxidoreductases
    peroxidase_keywords = ["hrp", "peroxidase"]
    # Metalloproteases
    mmp_keywords = ["mmp", "metalloproteinase"]
    # CD19 and B‑cell targets
    cd19_keywords = ["cd19", "b lymphocyte antigen"]

    # Helper to check if any keyword appears in the protein name
    def matches(keywords: List[str]) -> bool:
        return any(k in first_lower for k in keywords)

    if matches(gpcr_keywords):
        return (
            f"Research suggests that lipid composition strongly influences GPCR signalling; tailoring the phospholipid environment can stabilise {first} and enhance its activity."
        )
    if matches(ion_channel_keywords):
        return (
            f"Many ion channels are sensitive to their membrane environment; screening different nanodisc and detergent compositions can improve the stability and gating of {first}."
        )
    if matches(kinase_keywords):
        return (
            f"Kinase constructs often vary in solubility and activity; choosing appropriate domains and cofactors can dramatically increase yields and data quality for {first}."
        )
    if matches(tf_keywords):
        return (
            f"Transcription factors like {first} are typically disordered and benefit from solubility tags and specific additives to fold correctly."
        )
    if matches(disulfide_keywords):
        return (
            f"Disulfide‑rich proteins such as {first} require careful control of redox conditions and solubility tags to form correct disulfide bonds and remain active."
        )
    if matches(peroxidase_keywords):
        return (
            f"Peroxidases frequently require cofactors such as hemin, calcium and redox chaperones to fold properly; including these additives can dramatically increase yields."
        )
    if matches(mmp_keywords):
        return (
            f"Metalloproteinases often benefit from truncating non‑essential domains and supplementing zinc and calcium along with redox chaperones to improve expression."
        )
    if matches(cd19_keywords):
        return (
            f"Optimising the expression of {first} can involve screening full‑length and truncated constructs and using chaperones to assist folding."
        )
    # If we can't classify, ask a question to invite discussion
    return (
        f"I noticed your work on {first}.  Have you tried exploring different solubility tags or expression conditions to improve yields?"
    )


# ---------------------------------------------------------------------------
# Case study recommendation logic
#
# The following helper functions evaluate whether a case study from
# CASE_STUDIES_INFO is relevant to a lead based on the provided
# protein targets, therapeutic areas and topic.  A case study is
# considered relevant if the lead’s protein targets match any of the
# protein‑class keywords defined for the case, or if the therapeutic
# areas have any overlap with the case’s therapeutic‑area keywords.

def _tokenise(text: str) -> List[str]:
    """Split a comma‑ or space‑separated string into lower‑case tokens."""
    if not text:
        return []
    # Replace commas with spaces and split
    tokens = text.replace(",", " ").lower().split()
    return [t.strip() for t in tokens if t.strip()]


def recommend_case_studies(
    protein_targets: str,
    therapeutic_areas: str,
    topic: str | None = None,
) -> List[Dict[str, str]]:
    """
    Identify relevant case studies based on overlap with protein class or
    therapeutic area.

    Parameters
    ----------
    protein_targets : str
        Comma‑separated string of protein targets provided by the user.
    therapeutic_areas : str
        Comma‑separated string of therapeutic areas relevant to the lead’s work.
    topic : Optional[str]
        The mapped topic category (e.g. "membrane", "antibody", etc.).  If
        provided, the function may prioritise case studies aligned with this
        topic.  If None, all cases are considered.

    Returns
    -------
    List[Dict[str, str]]
        A list of dictionaries, each with keys ``name`` and ``summary`` for
        relevant case studies.  The order reflects the sequence in
        CASE_STUDIES_INFO.
    """
    # Flatten the input text into tokens for matching
    target_tokens = _tokenise(protein_targets)
    area_tokens = _tokenise(therapeutic_areas)
    # Prepare a combined string for simpler substring matching
    target_text = " ".join(target_tokens)
    area_text = " ".join(area_tokens)

    matches: List[Dict[str, str]] = []
    for case_name, info in CASE_STUDIES_INFO.items():
        # Optionally filter by topic: if a topic is provided, skip cases from
        # obviously unrelated categories.  For example, if the topic is
        # "membrane", skip soluble‑only examples like HRP unless therapeutic
        # area overlaps.  The mapping below defines which case names belong to
        # which topics.  Cases may appear in multiple topics (e.g. kinase is
        # relevant to both membrane and antibody due to shared workflows).
        case_topic: str | None = None
        # Determine a coarse topic for each case name
        if case_name in {"msba", "ffar1", "zmpste24", "bacteriorhodopsin", "dhodh"}:
            case_topic = "membrane"
        elif case_name in {"vegf", "hrp", "mmp1", "cd19", "tgfb1", "mdm2", "parkin", "vhl"}:
            case_topic = "generic"  # soluble / generic protein classes
        elif case_name in {"btk"}:
            case_topic = "kinase"
        elif case_name in {"gata"}:
            case_topic = "tf"
        # If a topic is given and doesn't match the case topic, only consider
        # the case when the protein class overlaps.  We avoid recommending
        # soluble cases for membrane leads (and vice versa) based solely on
        # broad therapeutic areas like "cancer".
        if topic and case_topic and topic != case_topic:
            # Defer cross‑topic recommendations unless there is a direct
            # protein class match
            overlap_class_cross = any(
                cls_kw in target_text for cls_kw in info["protein_class"]
            )
            if not overlap_class_cross:
                continue

        # Check for protein class overlap
        overlap_class = any(
            cls_kw in target_text for cls_kw in info["protein_class"]
        )
        overlap_area = any(
            area_kw in area_text for area_kw in info["therapeutic_areas"]
        )
        if overlap_class or overlap_area:
            # Look up a human‑readable title; fall back to the key if not found
            case_title = CASE_TITLES.get(case_name, case_name)
            matches.append({"name": case_title, "summary": info["summary"]})
    return matches


def generate_recommendation_report(
    title: str,
    company: str,
    topic_input: str,
    protein_targets: str = "",
    therapeutic_areas: str = "",
    linkedin_summary: str = "",
    first_name: str = "",
    last_name: str = "",
) -> Dict[str, object]:
    """
    Generate a report recommending case studies and resources based on lead
    information.

    Parameters
    ----------
    title : str
        The lead’s job title.
    company : str
        The lead’s company or organisation name.
    topic_input : str
        Description of the research area.
    protein_targets : str, optional
        Comma‑separated list of protein targets.
    therapeutic_areas : str, optional
        Comma‑separated list of therapeutic areas relevant to the lead’s work.
    linkedin_summary : str, optional
        Free‑text LinkedIn summary.

    Returns
    -------
    Dict[str, object]
        A dictionary containing the mapped topic, current customer flag, install
        location (if applicable), a list of recommended case studies with
        summaries, a list of suggested collateral items, and any prior
        campaign participation.  The dictionary also includes the raw
        segment and persona for potential downstream use.
    """
    # Infer segment and persona based on company, title and LinkedIn summary
    segment, persona = infer_segment_persona(company, title, linkedin_summary)

    # Determine workflow/topic.  If blank, default to generic.  Combine with summary to improve mapping.
    workflow_input = topic_input.strip() if topic_input else ""
    if workflow_input:
        topic = map_topic(f"{workflow_input} {linkedin_summary}")
    else:
        topic = "generic"

    # Override topic to antibody if therapeutic areas suggest biologics/antibodies
    areas_norm = _norm(therapeutic_areas)
    if areas_norm:
        for kw in ANTIBODY_KEYWORDS:
            if kw in areas_norm:
                topic = "antibody"
                break
        else:
            if "biologic" in areas_norm or "biologics" in areas_norm:
                topic = "antibody"

    # Identify if company is part of the installed base
    install_location = get_install_location(company)
    current_customer = install_location is not None

    # Find relevant case studies
    cases = recommend_case_studies(protein_targets, therapeutic_areas, topic)
    # Determine suggested collateral items
    resources = RESOURCE_SUGGESTIONS.get(topic, RESOURCE_SUGGESTIONS["generic"])

    # Look up campaigns associated with first/last name
    campaigns = []
    if first_name and last_name:
        campaigns = get_campaign_history(first_name, last_name)

    return {
        "segment": segment,
        "persona": persona,
        "topic": topic,
        "current_customer": current_customer,
        "install_location": install_location,
        "recommended_cases": cases,
        "suggested_resources": resources,
        "campaign_history": campaigns,
    }


# ---------------------------------------------------------------------------
# Email generation

# The email generation function has been removed.  The current focus of this
# module is to provide data and recommend appropriate case studies rather than
# generate outbound email sequences.  Previously this file contained a
# `generate_email_campaign` function that constructed personalised email
# sequences.  If email generation is needed in future, a new function can
# be added here to create outreach messages using the recommendation report.


# ---------------------------------------------------------------------------
# Command‑line interface

def main() -> None:
    """
    If executed as a script, prompt the user for input and display a
    recommendation report.  The report includes the mapped topic,
    segment/persona, a platform overview, a personalised research insight,
    relevant case studies and suggested collateral.  This function no
    longer generates email sequences; instead, it focuses on sharing
    data and recommending appropriate examples from Nuclera’s case study
    library.
    """
    import sys
    import textwrap

    print("Nuclera Case Study Recommendation Tool")
    print("Enter the lead’s information when prompted.  Leave a field blank if unknown.")

    try:
        first_name = input("First name: ").strip()
        last_name = input("Last name: ").strip()
        title = input("Lead title (e.g. Director of Protein Engineering): ").strip()
        company = input("Company/organisation: ").strip()
        workflow = input(
            "Workflow (soluble, membrane, antibody), leave blank for generic: "
        ).strip()
        targets = input("Protein targets (comma‑separated, optional): ").strip()
        areas = input("Therapeutic areas (comma‑separated, optional): ").strip()
        summary = input("LinkedIn summary (free text, optional): ").strip()
    except (EOFError, KeyboardInterrupt):
        print("\nAborted by user.")
        sys.exit(0)

    report = generate_recommendation_report(
        title=title,
        company=company,
        topic_input=workflow,
        protein_targets=targets,
        therapeutic_areas=areas,
        linkedin_summary=summary,
        first_name=first_name,
        last_name=last_name,
    )

    print("\nRecommendation report:\n")
    print("Workflow:", report["topic"])
    print("Segment / persona:", f"{report['segment']} / {report['persona']}")
    # Indicate whether this is a current customer
    if report["current_customer"]:
        loc = report.get("install_location", "") or "unknown location"
        print(f"\nCurrent customer: Yes (unit located in {loc})")
    else:
        print("\nCurrent customer: No")
    # List recommended case studies
    if report["recommended_cases"]:
        print("\nRecommended case studies:")
        for case in report["recommended_cases"]:
            print(f"  • {case['name']}: {case['summary']}")
    else:
        print("\nNo direct case study overlap found.  Consider sharing general workflow resources.")
    # List suggested collateral
    if report["suggested_resources"]:
        print("\nSuggested collateral:")
        for item in report["suggested_resources"]:
            print(f"  • {item}")
    # Show prior campaign participation if available
    if report["campaign_history"]:
        print("\nPrior campaign participation:")
        for camp in report["campaign_history"]:
            print(f"  • {camp['campaign_name']} (first associated on {camp['date']})")


if __name__ == "__main__":
    main()
