function dset = CollinsEtAl_datasets

% Copyright (c) 2012 The Regents of the University of California
% All Rights Reserved.

% This function sets a variable called dset, short for datasets, which
% contains category names, stimulus names, and other details about material
% used in seven tonal priming experiments. See the paper cited in the
% script CollinsEtAl_analysis for more details.

% Tom Collins, 2011.10.28.

nd = 0;
dset = struct();

nd = nd+1;
dset(nd).id = 'marmel_and_tillmann_mp_2009';
dset(nd).description = ['Marmel and Tillmann (2009) studied timbre and'...
    ' intonation tasks for mediant and leading-tone targets. File'...
    ' names containing iii refer to mediant stimuli, and vii refers' ...
	' to leading-tone stimuli. File names containing 35 could be'...
    ' ignored, as these are the mistuned stimuli (foils). There are 48'...
    ' stimuli in total.'];
dset(nd).stimCategoryLabels = {'Mediant', 'LeadingTone'};
dset(nd).stimCategoryLabelsShort = {'III', 'VII'};
dset(nd).stimCategoryMembership.Mediant = {'Riii0'};
dset(nd).stimCategoryMembership.LeadingTone = {'Rvii0'};
dset(nd).event_onsets = 0:789.5:5526.5;
% NB stimuli have different rhythmic profiles, so provided eighth-note
% onsets.

nd = nd+1;
dset(nd).id = 'marmel_etal_jep_2010';
dset(nd).description = ['Marmel, Tillmann and Delbe (2010) compare'...
    ' tonic and subdominant targets over two experiments, one with'...
    ' piano timbre and the other with pure tones. File names containing'...
    ' i refer to tonic stimuli, and iv refers to subdominant stimuli.'...
    ' There are 48 stimuli in total (24 piano and 24 pure tone). For'...
    ' the piano stimuli, there are bright and dull target timbres. For'...
    ' the pure-tone, there are two different types of target timbre, A'...
    ' and B.'];
dset(nd).stimCategoryLabels = {'PianoSubdominantBright',...
    'PianoTonicDull', 'PureTonicA', 'PureSubdominantA'};
dset(nd).stimCategoryLabelsShort = {'PnIVB', 'PnID', 'PuIA', 'PuIVA'};
dset(nd).stimCategoryMembership.PianoSubdominantBright = {'piano2BivF'};
dset(nd).stimCategoryMembership.PianoTonicDull = {'piano1MiF'};
dset(nd).stimCategoryMembership.PureTonicA = {'pure1AiF'};
dset(nd).stimCategoryMembership.PureSubdominantA = {'pure2AivF'};
dset(nd).event_onsets = 0:789.5:5526.5;
% NB stimuli have different rhythmic profiles, so provided eighth-note
% onsets.

nd = nd+1;
dset(nd).id = 'tillmann_etal_cbr_2003';
dset(nd).description = ['Tillmann, Janata and Bharucha (2003) studied'...
    ' a dissonance discrimination task for related and unrelated'...
    ' targets. As an example, the filename 1RCinBB should be read as'...
    ' Related target, Consonant, in B major. The final B is just a'...
    ' label. Details: Sixteen eight-chord sequences were selected from'...
    ' Bigand and Pineau (1997), Pineau and Bigand (1997), and Bigand,'...
    ' Tillmann, Poulin, and D`Adamo (2001) that ended on the tonic.'...
    ' Duplicates were created in C and B major. Duplicates in C major'...
    ' were also created where the final chord was transposed down a'...
    ' semitone (giving B major), and duplicates were created in B major'...
    ' where the final chord was transposed up a semitone (giving C'...
    ' major). Thus the related target chord in C major matched the'...
    ' unrelated target chord in B major, and vice versa. Finally,'...
    ' consonant and dissonant versions of each sequence were created,'...
    ' giving 128 stimuli overall.'];
dset(nd).stimCategoryLabels = {'RelatedConsonantInB','RelatedConsonantInC', ...
    'UnrelatedConsonantInB', 'UnrelatedConsonantInC'};
dset(nd).stimCategoryLabelsShort = {'Rel B', 'Rel C','Un B', 'Un C'};
dset(nd).stimCategoryMembership.RelatedConsonantInB = {'1RCinBB'};
dset(nd).stimCategoryMembership.RelatedConsonantInC = {'1RCinCB'};
dset(nd).stimCategoryMembership.UnrelatedConsonantInB = {'3UCinBB'};
dset(nd).stimCategoryMembership.UnrelatedConsonantInC = {'3UCinCB'};
dset(nd).event_onsets = 0:500:3500;

nd = nd+1;
dset(nd).id = 'tillmann_etal_jep_2008';
dset(nd).description = ['Tillmann, Janata, Birk and Bharucha (2008)'...
    ' studied a dissonance discrimination task for tonic (T), dominant'...
    ' (D), subdominant (S), and baseline (b) targets. As an example,'...
    ' the filename ADb should be read as a stimulus ending on the'...
    ' dominant, with a baseline context. The initial A is just a label.'...
    ' Details: These stimuli were split algorithmically from a single'...
    ' wav file by Charles Delbé. Originally they were twelve eight-'...
    ' chord sequences selected from Bigand and Pineau (1997) and Pineau'...
    ' and Bigand (1997) that ended on the tonic. The penultimate chord'...
    ' was altered to G2, G3, D4, and D5. All sequences in this set end'...
    ' this way. Different contexts were then composed so that the'...
    ' penultimate and target chords evoked ii-V progressions (labelled'...
    ' dominant) and I-IV progressions (labelled subdominant). For each'...
    ' of the twelve stimuli from the three categories (tonic, dominant,'...
    ' subdominant), a baseline sequence was composed that crisscrossed'...
    ' the cycle of fifths, giving 72 stimuli in total.'];
dset(nd).stimCategoryLabels = {'Tonic', 'Dominant', 'Subdominant',...
    'TonicBaseline', 'DominantBaseline', 'SubdominantBaseline'};
dset(nd).stimCategoryLabelsShort = {'I', 'V', 'IV',...
    'I BL', 'V BL', 'IV BL'};
dset(nd).stimCategoryMembership.Tonic = {'CT'};
dset(nd).stimCategoryMembership.Dominant = {'CD'};
dset(nd).stimCategoryMembership.Subdominant = {'CS'};
dset(nd).stimCategoryMembership.TonicBaseline = {'CTb'};
dset(nd).stimCategoryMembership.DominantBaseline = {'CDb'};
dset(nd).stimCategoryMembership.SubdominantBaseline = {'CSb'};
dset(nd).event_onsets = 0:500:3500;

nd = nd+1;
dset(nd).id = 'tillmann_etal_jep_2003';
dset(nd).description = ['Tillmann, Janata, Birk and Bharucha (2003)'...
    ' studied a dissonance discrimination task for tonic (T),'...
    ' subdominant (S), and baseline (b) targets (experiment 2). As an'...
    ' example, the filename J_IV_baseline_cons should be read as'...
    ' stimulus ending consonantly on the subdominant with a baseline'...
    ' context. The initial J is just a label. Details: There appear to'...
    ' be some extra stimuli here, such as I__L_cons.wav. Originally'...
    ' they were twelve eight-chord sequences selected from Bigand and'...
    ' Pineau (1997) and Pineau and Bigand (1997) that ended on the'...
    ' tonic, and twelve that ended on the subdominant. All twelve major'...
    ' keys are represented within each category. For each of the twelve'...
    ' stimuli from the two categories (tonic, subdominant), a baseline'...
    ' sequence was composed that crisscrossed the cycle of fifths,'...
    ' giving 48 stimuli. The other 48 stimuli are dissonant foils, so'...
    ' could be ignored. These stimuli were provided by Barbara Tillmann'...
    ' on 22 December 2011 and were split into individual stimuli by Tom'...
    ' Collins using the function mironsets (see splitting_stimuli.m)'];
dset(nd).stimCategoryLabels = {'TonicConsonant',...
    'TonicBaselineConsonant', 'SubdominantConsonant',...
    'SubdominantBaselineConsonant'};
dset(nd).stimCategoryLabelsShort = {'I', 'I BL', 'IV', 'IV BL'};
dset(nd).stimCategoryMembership.TonicConsonant = {'C_I_cons'};
dset(nd).stimCategoryMembership.TonicBaselineConsonant =...
    {'C_I_baseline_cons'};
dset(nd).stimCategoryMembership.SubdominantConsonant = {'C_IV_cons'};
dset(nd).stimCategoryMembership.SubdominantBaselineConsonant =...
    {'C_IV_baseline_cons'};
dset(nd).event_onsets = 0:500:3500;

nd = nd+1;
dset(nd).id = 'bigand_etal_jep_2003';
dset(nd).description = ['Bigand, Poulin, Tillmann, Madurell, and'...
    ' D`Adamo (2003) studied a dissonance discrimination task for tonic'...
    ' (I) and subdominant (IV) targets, some of whose contexts'...
    ' contained the target chord and some of whose contexts contained'...
    ' the subdominant chord. As an example, the filename K_I_NTC should'...
    ' be read as a stimulus ending on the tonic with No Target in'...
    ' Context, whereas K_I_CIV should be read as stimulus ending on the'...
    ' tonic with Context containing IV. The initial K is just a label.'...
    ' Details: These stimuli were split algorithmically by Tom Collins'...
    ' from a single wav file supplied by Charles Delbé. There are'...
    ' twelve eight-chord sequences that end on the tonic but whose'...
    ' context does not contain the tonic. These twelve are manipulated'...
    ' to create another twelve stimuli that end on the tonic and whose'...
    ' contexts contain one or two occurrences of the subdominant. For'...
    ' each of the previously mentioned twelve stimuli, a version is'...
    ' created that ends on the subdominant. And for each of these'...
    ' twelve stimuli, a version is created where the context does not'...
    ' contain the target (subdominant), giving 48 stimuli in total.'];
dset(nd).stimCategoryLabels = {'TonicNTC', 'SubdominantNTC'};
dset(nd).stimCategoryLabelsShort = {'I NTC', 'I CIV', 'IV NTC', 'IV CIV'};
dset(nd).stimCategoryMembership.TonicNTC = {'B_I_NTC'};
dset(nd).stimCategoryMembership.SubdominantNTC = {'B_IV_NTC'};
dset(nd).event_onsets = 0:600:4200;

nd = nd+1;
dset(nd).id = 'marmel_etal_pp_2008';
dset(nd).description = ['Marmel et al. (2008) studied intonation tasks'...
    ' for melodic tonic and subdominant targets. File names containing'...
	' i refer to tonic stimuli, and iv refers to subdominant stimuli.'...
    ' File names containing 35 could be ignored, as these are the'...
    ' mistuned stimuli (foils). There are 48 stimuli in total.'];
dset(nd).stimCategoryLabels = {'Tonic', 'Subdominant'};
dset(nd).stimCategoryLabelsShort = {'I', 'IV', 'I mis', 'IV mis'};
dset(nd).stimCategoryMembership.Tonic = {'Gi'};
dset(nd).stimCategoryMembership.Subdominant = {'Giv'};
dset(nd).event_onsets = 0:789.5:5526.5;
% NB stimuli have different rhythmic profiles, so provided eighth-note
% onsets.

end