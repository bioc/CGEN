%macro table2Format(data=_LAST_, idVar=, outfile="out", delimiter="|",
                    outMiss=-9999, startCol=1, stopCol=-1, idFile="ids.txt",
                    vars=);
  /* This macro will first change the missing values for all numeric variables
	 to outMiss. Then it will convert the character variables into numeric
	 variables, again changing the missing values to outMiss. The new codes for
	 the character variables will be 1, 2, 3, ....
	 The data will then be transposed and exported.

     History:  Apr 15 2008 Let outMiss be the missing value for both
					       numeric and character variables.


    data         The sas data set containing the snp variables and the 
                 subject id variable.
                 The default is _LAST_
    idVar        The id variable. No default.
    outfile      The output file. Must be quoted.
                 The default is "out".
    delimiter    The delimiter in the output file. Must be quoted.
                 The default is "|".
	outMiss      Numeric value to denote missing values.
			     The default is -9999.
    startCol     Starting column to use.
			     The default is 1
    stopCol      Stopping column. The default is -1.
    idFile       Name of the file that will contain the ids.
				 The default is "id.txt".
	vars         Variables to keep. This overrides start.col and stop.col
  */

  %local temp nv tempData i idFlag error j missChar useData useFlag transData;
  %local idData dropStr nextData varFlag n_num n_char nv2 nobs k;

  /* Initialize */
  %let tempData  = _tmp_tab2form_4134;
  %let idFlag    = %length(&idVar);
  %let error     = 0;
  %let useData   = _use_tab2form_4134;
  %let idData    = _id_tab2form_4134;
  %let useFlag   = 0;
  %let varFlag   = %length(&vars);
  %let nextData  = &data;
  %let n_num     = 0;
  %let n_char    = 0;
  %if (%length(&outMiss) = 0) %then %let outMiss = -9999;
  %let missChar  = "&outMiss";

  %if (&idFlag) %then %do;
    /* Check the id variable */
    %checkVar(&idVar, data=&data, _err=error, errMsg=0);
    %if (&error) %then %goto exit;

    %let dropStr = (drop=&idVar);
  %end;
  %else %let dropStr = ;

  /* Get the ids */
  data &idData; set &data(keep=&idVar); run;

  /* Get the number of variables */
  %getAttrNC(&data, inNList=nvars, _outNList=nv);

  /* Check stop.col */
  %if (&stopCol < 1) %then %let stopCol = &nv; 

  /* See if variables were specified */
  %if (&varFlag) %then %do;
    data &useData; set &data&dropStr;
	  keep &vars;
	run;

    %let nextData = &useData;
	%let startCol = 1;
	%let stopCol  = &nv;

	%let error = 0;
    %checkVar(&idVar, data=&useData, _err=error, errMsg=0);
    %if (&error) %then %let dropStr = ;
    %else              %let dropStr = (drop=&idVar);
  %end;

  %if (&startCol > 1) | (&stopCol < &nv) %then %do;
    /* Get the variable names. */
    proc contents data=&data&dropStr noprint
      out=&tempData(keep=NAME VARNUM);
    run;

    /* Sort by variable number */
	proc sort data=&tempData;
	  by VARNUM;
	run;

    /* Get the subset of variable names */
	%let temp = %eval(&stopCol - &startCol + 1);
	%let temp = %eval(&startCol + &temp - 1);
	data &tempData; 
      set &tempData(firstobs=&startCol obs=&temp);
	run;

    /* Get the number of observations (variables) */
    %getAttrNC(&tempData, inNList=nobs, _outNList=nv);

    /* Define local macro variables */
    %do i = 1 %to &nv;
      %local var&i; 
    %end;

	/* Save variable names in macro variables */
    %getMacroVars(&tempData, NAME, _mPrefix=var); 

    /* Subset the data */
    data &useData; set &data&dropStr;
      keep
      %do i = 1 %to &nv;
        &&var&i
	  %end;
	  ;
	run;
	%let dropStr  = ;
	%let useFlag  = 1;
	%let nextData = &useData

    /* Free macro variables */
    %freeMacroVars(list=, prefixList=var, lenList=&nv); 
  %end;

  /* Get the variable names. Only numeric variables */
  proc contents data=&nextData&dropStr noprint
    out=&tempData(keep=NAME TYPE where=(TYPE=1));
  run;

  /* Get the number of observations (variables) */
  %getAttrNC(&tempData, inNList=nobs, _outNList=nv);

  %if (&nv) %then %do;
    %let n_num = &nv;

    /* Define local macro variables */
    %do i = 1 %to &nv;
      %local var&i; 
    %end;

	/* Save variable names in macro variables */
    %getMacroVars(&tempData, NAME, _mPrefix=var); 

    /* Change missing values */
    data &useData; set &nextData;
	  %do i = 1 %to &nv;
        if (missing(&&var&i)) then &&var&i = &outMiss;
	  %end;
    run;

    /* Free macro variables */
    %freeMacroVars(list=, prefixList=var, lenList=&nv);

	%let useFlag  = 1;
	%let nextData = &useData;
  %end;

  /* Get the variable names. Only character variables */
  proc contents data=&nextData&dropStr noprint
    out=&tempData(keep=NAME TYPE VARNUM where=(TYPE=2));
  run;

  /* Get the number of observations */
  %getAttrNC(&tempData, inNList=nobs, _outNList=nv);
  %let n_char = &nv;

  /* Make sure that the variables are either all numeric or all character */
  %if ((&n_num) & (&n_char)) %then %do;
    %put ERROR: The variables must either be all numeric or all character;
	%let error = 1;
	%goto exit;
  %end;

  /* If no character variables, go straight to transposing */
  %if (&nv = 0) %then %do;
    %let transData = &nextData;
    %goto transpose;
  %end;

  %let transData = &tempData;

  /* Sort by variable number */
  proc sort data = &tempData;
    by VARNUM;
  run;

  /* Define local macro variables */
  %do i = 1 %to &nv;
    %local var&i old&i new&i n&i; 
  %end;

  /* Save variable names in macro variables */
  %getMacroVars(&tempData, NAME, _mPrefix=var);

  /* Get the subset of the data, and change missing values */
  %let error = 0;
  data &tempData; set &nextData;
    keep 
    %do i = 1 %to &nv;
      &&var&i
    %end;
	;
	%do i = 1 %to &nv;
	  if (&&var&i = &missChar) then do;
        call symput("error", 1);
		stop;
	  end;
      if (missing(&&var&i)) then &&var&i = &missChar;
	%end;
  run;

  /* Check for error */
  %if (&error) %then %do;
    %put ERROR: (Internal error) Define a new character for missing values!!!;
	%goto exit;
  %end;

  /* Call proc means */
  proc means data=&tempData noprint;
    class _all_;
    types
	%do i = 1 %to &nv;
      &&var&i 
	%end;
    ;
    output out = &tempData;
  run;

  /* Get the codes for the new numeric variables */
  data _NULL_; 
    length _var $32;
    set &tempData;
	retain _i 1 _var " ";
    %do i = 1 %to &nv;
      if (missing(&&var&i) = 0) then do;
	    if (_var ^= "&&var&i") then do;
          _i   = 1;
		  _var = "&&var&i";
		end;
        call symput(trim("old&i._") || left(_i), trim(left(&&var&i)));
        call symput(trim("new&i._") || left(_i), _i);
        call symput("n&i", trim(left(_i)));
        _i = _i + 1;
	  end;
    %end; 
  run;

  /* Create the new data set */
  %let error = 0;
  data &tempData; set &nextData;
    drop
	%do i = 1 %to &nv;
      &&var&i
	%end;
    ;
  
    %do i = 1 %to &nv;
	  %let temp = &&var&i;
	  if (missing(&temp)) then &temp._ = &outMiss;

	  %do j = 1 %to &&n&i;
        else if (&temp = "&&old&i._&j") then &temp._ = &&new&i._&j;
	  %end;
	  else do;
		call symput("error", 1);
		stop;
	  end;
	%end;
  run;

  /* Check for error */
  %if (&error = 1) %then %do;
    %put ERROR: in assigning numeric values;
	%goto exit;
  %end;

  /* Rename the variables */
  proc datasets library=work nodetails;
    modify &tempData;
	rename
	%do i = 1 %to &nv;
	  %let temp = &&var&i;
      &temp._ = &temp 
	%end;
	;
  run;
  quit;

  /* Free macro variables */
  %freeMacroVars(list=, prefixList=var, lenList=&nv);

  %transpose:

  /* Check the id variable */
  %let error = 0;
  %checkVar(&idVar, data=&transData, _err=error, errMsg=0);
  %if (&error) %then %let dropStr = ;
  %else              %let dropStr = (drop=&idVar);

  /* Transpose the data */
  proc transpose data=&transData&dropStr out=&tempData;
  run;
 
  /* Change the character variables back */
  %if (&n_char) %then %do;
    /* Get the number of observations */
    %getAttrNC(&tempData, inNList=nobs nvars, _outNList=nobs nv2);
	%let nv2 = %eval(&nv2 - 1);

    %let error = 0;
    data &tempData; set &tempData;
      drop col1 - col&nv2;
      select (_N_);
	  %do k = 1 %to &nobs;
	    when (&k) do;
        %do i = 1 %to &nv2;
	      if (col&i = &missChar) then col&i._ = &missChar;
	      %do j = 1 %to &&n&k;
            else if (col&i = &&new&k._&j) then col&i._ = "&&old&k._&j";
	      %end;
	      else do;
		    call symput("error", 1);
		    stop;
	      end;
		 %end;
        end;
	  %end;
	  end;
    run;

    /* Check for error */
    %if (&error) %then %do;
      %put ERROR: Changing back to characters;
	  %goto exit;
    %end;

  %end; /* END: %if (&n_char) %then %do  */ 

  /* Free macro variables */
  %do i = 1 %to &nv;
    %do j = 1 %to &&n&i;
      %let n&i = ;
	  %let old&i._&j = ;
	  %let new&i._&j = ;
	%end;
  %end;

  %let error = 0;
  %checkVar(_LABEL_, data=&tempData, _err=error, errMsg=0);
  %if (&error) %then %let temp = ;
  %else              %let temp = (drop=_LABEL_);

  /* Export the data */
  %export(data=&tempData&temp, outfile=&outfile, delimiter=&delimiter, header=0);
  %export(data=&idData, outfile=&idFile, delimiter=" ", header=0);

  %exit: ;

  %deleteData(&tempData);
  %deleteData(&idData);
  %if (&useFlag) %then %deleteData(&useData);

%mend table2Format;

%macro export(data=_LAST_, outfile="out.txt", delimiter='09'x, header=1);
  /* Macro to export a sas data set to a file
     data        Sas data set. The default is _LAST_
     outfile     Output file (in quotations). The default is "out.txt"
     delimiter   For tab delimited files, use '09'x, otherwise put the
                 delimiter in quotations.
                 The default is '09'x
     header      0 or 1 to include the variable names in the exported file.
                 The default is 1.
  */
 
  %local nv tempData i temp nvm1;

  /* Initialize */
  %let tempData  = _tmp_export_9754;

  /* Get the variable names */
  proc contents data=&data noprint
    out=&tempData(keep=NAME TYPE LENGTH VARNUM);
  run;
  proc sort data = &tempData;
    by VARNUM;
  run;

  /* Get the number of observations */
  %getAttrNC(&tempData, inNList=nobs, _outNList=nv);

  /* Define local macro variables */
  %do i = 1 %to &nv;
    %local var&i type&i len&i; 
  %end;

  /* Save variable names in macro variables */
  %getMacroVars(&tempData, NAME TYPE LENGTH, _mPrefix=var type len);

  /* Change the macro variables for the data step */
  %do i = 1 %to &nv;
    %if (&&type&i = 1) %then %do;
       %let type&i = best12.;
       %let len&i  = ;
	%end;
	%else %do;
      %let type&i = $&&len&i...;
	  %let len&i = $;
	%end;
  %end;

  %let nvm1 = %eval(&nv - 1);

  /* Export */
  data _NULL_; set &data;
    file &outfile DSD STOPOVER lrecl=32767 delimiter=&delimiter;     
	/* Set the formats */
    %do i = 1 %to &nv;
      format &&var&i &&type&i;
	%end;
	/* Print the header */
	%if (&header) %then %do;
      if (_N_ = 1) then do;
	    put
        %do i = 1 %to &nvm1;
          "&&var&i"
          &delimiter
	    %end;
        "&&var&nv"
	    ;
	  end;
	%end;
	/* Print the data */
	do;
      %do i = 1 %to &nvm1;
        put &&var&i &&len&i @;
	  %end;
	  put &&var&nv &&len&i;
	  ;
	end;

  run;

  %exit: ;

  %deleteData(&tempData);

%mend export;

%macro checkDataName(data, _err=error, exist=0);
  /* Macro to check the validity of a sas data set name.
     Returns 1 if an error was encountered 
     data   Sas data set to check.
     _err   Name of the macro error variable.
     exist  0 or 1 for determining if the data set exists.
  */

  %local __lib __file __temp;

  /* Get the file name */
  %getFilename(&data, _filename=__file);
  %let __temp = %sysfunc(nvalid(&__file));
  %if (&__temp = 0) %then %do;
    %put ERROR: The data set name &data is invalid;
    %let &_err = 1;
  %end;

  /* Get the library name */
  %getLibname(&data, _libname=__lib);

  %if (%sysfunc(libref(&__lib))) %then %do;
    %put ERROR: For the data set &data, the library &__lib is not defined;
    %let &_err = 1;
  %end;

  %if (&exist) %then %do;
    %if (%sysfunc(exist(&data)) = 0) %then %do;
      %put ERROR: The data set &data does not exist;
      %let &_err = 1;
    %end;
  %end;

%mend checkDataName;

%macro checkVar(var, dsid=0, data=, _err=error, vartype=, _retVarlen=,
                errMsg=1);
    /* Macro to check the existence of a variable on a sas data set, or checks
       the validity of a sas variable name.
       dsid        File id from opening the data set.
                   Use dsid=-1, to only check for a valid variable name.
       var         Variable to check.
       data        Name of the data set.
       _err        Name of macro error variable.
       vartype     Type of variable (N or C). If not specified, then the
                   variable type will not be checked. 
       _retVarlen  Returned variable length
	   errMsg      0 or 1 to print an error message
    */
    %local __temp __flag ;

    /* Check variable name */
	%if (&dsid = -1) %then %do;
      %let __temp = %sysfunc(nvalid(&var));
      %if (&__temp = 0) %then %do;
        %if (&errMsg) %then %put ERROR: The variable name &var is invalid;
        %let &_err = 1;
      %end;
	  %goto exit_checkVar;
	%end;

    /* Get the data set id if not passed in */
	%if (&dsid = 0) %then %do;
      %let dsid   = %sysfunc(open(&data));
	  %let __flag = 1;
	%end;
	%else %let __flag = 0;

    /* The varnum function is used to return the column number of a variable */
    %let __temp = %sysfunc(varnum(&dsid, &var));
    %if (&__temp <= 0) %then %do;
	  %let &_err = 1;
      %if (&errMsg) %then %put ERROR: The variable &var was not found in the data set &data;
    %end;

    /* Get the variable length */
    %if (%length(&_retVarlen)) %then %do;
      %let &_retVarlen = %sysfunc(varlen(&dsid, &__temp)); 
	%end;

	/* Check the variable type */
    %if (%length(&vartype)) %then %do;
	  %let __temp  = %sysfunc(vartype(&dsid, &__temp));
	  %let __temp  = %upcase(%substr(&__temp, 1, 1));
	  %let vartype = %upcase(%substr(&vartype, 1, 1));
	  %if (&__temp ^= &vartype) %then %do;
        %let &_err = 1;
        %if (&errMsg) %then %put ERROR: The variable &var is not the proper type;
	  %end;
    %end;

    %if (&__flag) %then %let __temp = %sysfunc(close(&dsid));

    %exit_checkVar: 

%mend checkVar;

%macro deleteData(data, type=DATA);
  /* Macro to delete a data set.
     data   Sas data set to be deleted. No default.
     type   DATA or VIEW
  */
  %if %sysfunc(exist(&data, &type)) %then %do;
    %if (%upcase(&type) ^= VIEW) %then %let type = TABLE;
    proc sql; drop &type &data; quit;
  %end;
%mend deleteData;

%macro freeMacroVars(list=, prefixList=, lenList=);
  /* Macro to free macro variables.
     list         List of macro variables seperated by a space.
     prefixList   List of macro prefixes
     lenList      Number of macros for each prefix.
     Example: freeMacroVars(list=a b c, prefixList=x y, lenList=3 2) will free the
     macro variables a, b, c, x1, x2, x3, y1, y2
  */
  %local __stop __i __temp __temp2 __j;

  %if (%length(&list)) %then %do;
    %let __stop = 0;
    %let __i    = 1;
    %do %until (&__stop = 1);
      /* Use the scan function to pick off each element */
      %let __temp = %scan(&list, &__i, ' ');
      %if %length(&__temp) > 0 %then %do;
        %let &__temp = ;
        %let __i = %eval(&__i + 1);
	  %end;
	  %else %let __stop = 1;
    %end;
  %end;

  %if ((%length(&prefixList)) and (%length(&lenList))) %then %do;
    %let __stop = 0;
    %let __i    = 1;
    %do %until (&__stop = 1);
      /* Use the scan function to pick off each element */
      %let __temp  = %scan(&prefixList, &__i, ' ');
      %let __temp2 = %scan(&lenList, &__i, ' ');
      %if ((%length(&__temp)) and (%length(&__temp2))) %then %do;
        %do __j = 1 %to &__temp2;
          %let &__temp&__j = ;
		%end;
        %let __i = %eval(&__i + 1);
	  %end;
	  %else %let __stop = 1;
    %end;
  %end;

%mend freeMacroVars;

%macro getAttrNC(data, dsid=0, inNList=, inCList=, _outNList=, _outCList= );
  /* Macro to get data set info.
     data         Sas data set
     dsid         Data set id (if exists)
     inNList      List of numeric attributes
     inCList      List of character attributes
     _outNList    Output list of macro variable names for the numeric attributes
     _outCList    Output list of macro variable names for the character attributes
  */
  %local __rc __flag __stop __tempIn __i __tempOut;

  %if (&dsid = 0) %then %do;
    %let dsid   = %sysfunc(open(&data));
	%let __flag = 1;
  %end;
  %else %let __flag = 0;

  /* For numeric items */
  %if ((%length(&inNList) > 0) and (%length(&_outNList) > 0)) %then %do;
    %let __stop = 0;
    %let __i    = 1;
    %do %until(&__stop = 1);
      %let __tempIn  = %scan(&inNList, &__i, ' ');
      %let __tempOut = %scan(&_outNList, &__i, ' ');
      %if ((%length(&__tempIn) > 0) and (%length(&__tempOut) > 0)) %then %do;
        %let &__tempOut = %sysfunc(attrn(&dsid, &__tempIn));
        %let __i        = %eval(&__i + 1);
  	  %end;
	  %else %let __stop = 1;
    %end;
  %end;

  /* For character items */
  %if ((%length(&inCList) > 0) and (%length(&_outCList) > 0)) %then %do;
    %let __stop = 0;
    %let __i    = 1;
    %do %until(&__stop = 1);
      %let __tempIn  = %scan(&inCList, &__i, ' ');
      %let __tempOut = %scan(&_outCList, &__i, ' ');
      %if ((%length(&__tempIn) > 0) and (%length(&__tempOut) > 0)) %then %do;
        %let &__tempOut = %sysfunc(attrc(&dsid, &__tempIn));
        %let __i        = %eval(&__i + 1);
	  %end;
	  %else %let __stop = 1;
    %end;
  %end;
  
  %if (&__flag) %then %let __rc = %sysfunc(close(&dsid));

%mend getAttrNC;

%macro getMacroVars(data, vars, _mPrefix=v_);
  /* Macro to set macro variables from a sas variable.
     data       Sas data set.
     vars       List of variables seperated by a blank space
     _mPrefix   List of macro variable prefixes that will be defined as &_mPrefix1,...,&_mPrefixn
                for i = 1,...,n observations
  */
  %local __stop __i __tempv __tempm;
  %let __stop = 0;
  %let __i    = 1;

  data _NULL_; set &data;
    %do %until(&__stop = 1);
      %let __tempv = %scan(&vars, &__i, ' ');
      %let __tempm = %scan(&_mPrefix, &__i, ' ');
      %if ((%length(&__tempv) > 0) and (%length(&__tempm) > 0)) %then %do;
        call symput(trim("&__tempm") || left(_N_), trim(left(&__tempv)));
		%let __i = %eval(&__i + 1);
	  %end;
	  %else %let __stop = 1;
	%end;
  run;

%mend getMacroVars;

