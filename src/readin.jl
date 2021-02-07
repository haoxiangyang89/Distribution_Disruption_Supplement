# functions to read-in the data

function readMP(fileAdd::String)
    # read in the file from a specific address
    # input: the address of the Matpower data file
    # output:
    #   busST: the string containing bus information
    #   genST: the string containing generator information
    #   brST: the string containing branch information
    #   cST: the string containing cost information
    #   baseMVA: the base unit

    f = open(fileAdd);
    rawStr = read(f,String);
    close(f);

    # separate the raw strings for bus/gen/branch/cost
    if typeof(match(r"mpc.bus = \[([0-9\t\ \.\;\n\r\-e]*)\];",rawStr)) != Nothing
        busST = match(r"mpc.bus = \[([0-9\t\ \.\;\n\r\-e]*)\];",rawStr).captures[1];
    else
        busST = "";
    end
    if typeof(match(r"mpc.gen = \[([0-9\t\ \.\;\n\r\-e]*)\];",rawStr)) != Nothing
        genST = match(r"mpc.gen = \[([0-9\t\ \.\;\n\r\-e]*)\];",rawStr).captures[1];
    else
        genST = "";
    end
    if typeof(match(r"mpc.branch = \[([0-9\t\ \.\;\n\r\-e]*)\];",rawStr)) != Nothing
        brST = match(r"mpc.branch = \[([0-9\t\ \.\;\n\r\-e]*)\];",rawStr).captures[1];
    else
        brST = "";
    end
    if typeof(match(r"mpc.gencost = \[([0-9\t\ \.\;\r\n\-e]*)\];",rawStr)) != Nothing
        cST = match(r"mpc.gencost = \[([0-9\t\ \.\;\n\r\-e]*)\];",rawStr).captures[1];
    else
        cST = "";
    end
    if typeof(match(r"mpc.uncertain = \[([0-9\t\ \.\;\n\r\-e]*)\];",rawStr)) != Nothing
        uST = match(r"mpc.uncertain = \[([0-9\t\ \.\;\n\r\-e]*)\];",rawStr).captures[1];
    else
        uST = ""
    end

    baseMVA = parse(Float64,match(r"mpc.baseMVA = ([0-9\.]+);",rawStr).captures[1]);

    return busST,genST,brST,cST,uST,baseMVA;
end

function parsebusST(busST,baseMVA)
    # parse the busST string collected from readMP
    # input: the strings containing bus information
    # output: List of bus ID, list of generator ID, data variable Vmax/Vmin

    # separate the raw strings line by line
    busSTlbl = collect(eachmatch(r"([0-9\t\ \.\-e]+[\;\r]*\n)",busST));

    # parse each line of the bus data
    Vmax = Dict();
    Vmin = Dict();
    bType = Dict();
    gs = Dict();
    bs = Dict();
    Vmag = Dict();
    Vang = Dict();
    Pd = Dict();
    Qd = Dict();
    IDList = [];
    for bstr in busSTlbl
        bsstr = bstr.captures[1];
        bsstr = strip(strip(bsstr),';');
        bdata = split(bsstr);

        # update the IDList with the ID of the current bus
        cID = parse(Int64,bdata[1]);
        push!(IDList,cID);

        # obtain the type of the bus
        bType[cID] = parse(Int64,bdata[2]);

        # obtain the active/reactive power demand
        Pd[cID] = parse(Float64,bdata[3])/baseMVA;
        Qd[cID] = parse(Float64,bdata[4])/baseMVA;

        # shunt conductance and susceptance
        gs[cID] = parse(Float64,bdata[5])/baseMVA;
        bs[cID] = parse(Float64,bdata[6])/baseMVA;

        # nominal voltage magnitude and angle
        Vmag[cID] = parse(Float64,bdata[8]);
        Vang[cID] = parse(Float64,bdata[9]);

        # obtain the maximum and the minimum of the voltage magnitude
        Vmax[cID] = parse(Float64,bdata[12]);
        Vmin[cID] = parse(Float64,bdata[13]);
    end

    return IDList,Vmax,Vmin,gs,bs,Vmag,Vang,Pd,Qd,bType;
end

function parsegenST(genST,baseMVA)
  # parse the genST string collected from readMP
  # input: the strings containing generator information
  # output: data variable

  # separate the raw strings line by line
    genSTlbl = collect(eachmatch(r"([0-9\t\ \.\-e]+[\;\r]*\n)",genST));

    # parse each line of the generator data
    Pmax = Dict();
    Pmin = Dict();
    Qmax = Dict();
    Qmin = Dict();
    Pg = Dict();
    Qg = Dict();
    RU = Dict();
    RD = Dict();
    L = Dict();
    LR = Dict();
    l = 0;
    genIDList = [];
    for gstr in genSTlbl
        gsstr = gstr.captures[1];
        gsstr = strip(strip(gsstr),';');
        gdata = split(gsstr);
        l += 1;

        cID = parse(Int64,gdata[1]);
        push!(genIDList,l)

        # obtain the maximum/minimum active/reactive power output
        # L is the bus in which the generator is located
        L[l] = cID;
        if !(cID in keys(LR))
            LR[cID] = [];
            push!(LR[cID],l);
        else
            push!(LR[cID],l);
        end
        Pmax[l] = parse(Float64,gdata[9])/baseMVA;
        Pmin[l] = parse(Float64,gdata[10])/baseMVA;
        Qmax[l] = parse(Float64,gdata[4])/baseMVA;
        Qmin[l] = parse(Float64,gdata[5])/baseMVA;

        Pg[l] = parse(Float64,gdata[2])/baseMVA;
        Qg[l] = parse(Float64,gdata[3])/baseMVA;

        RU[l] = parse(Float64,gdata[17]);
        RD[l] = -RU[l];
    end
    return genIDList,L,LR,Pmax,Pmin,Qmax,Qmin,Pg,Qg,RU,RD;
end

function parsebrST(brST,baseMVA)
    # parse the brST string collected from readMP
    # input: the strings containing branch information
    # output: data variable

    # separate the raw strings line by line
    brSTlbl = collect(eachmatch(r"([0-9\t\ \.\-e]+[\;\r]*\n)",brST));

    g = Dict();
    b = Dict();
    bc = Dict();
    angmax = Dict();
    angmin = Dict();
    τ1 = Dict();
    τ2 = Dict();
    σ = Dict();
    m = Dict();
    brList = [];
    brList1 = [];
    rateA = Dict();

    for brstr in brSTlbl
        brsstr = brstr.captures[1];
        brsstr = strip(strip(brsstr),';');
        brdata = split(brsstr);

        fID = parse(Int64,brdata[1]);
        tID = parse(Int64,brdata[2]);
        if !((fID,tID) in brList)
            m[(fID,tID)] = 1;
            m[(tID,fID)] = 1;
            push!(brList,(fID,tID));
            push!(brList,(tID,fID));

            r = parse(Float64,brdata[3]);
            x = parse(Float64,brdata[4]);
            bc1 = parse(Float64,brdata[5]);
            τ = parse(Float64,brdata[9]);
            if τ == 0.0
                τ = 1.0;
            end
            σ1 = parse(Float64,brdata[10]);
            g[(fID,tID,1)] = round(r/(r^2 + x^2),digits=6);
            g[(tID,fID,1)] = round(r/(r^2 + x^2),digits=6);
            b[(fID,tID,1)] = round(-x/(r^2 + x^2),digits=6);
            b[(tID,fID,1)] = round(-x/(r^2 + x^2),digits=6);
            bc[(fID,tID,1)] = bc1;
            bc[(tID,fID,1)] = bc1;
            τ1[(fID,tID,1)] = τ;
            τ1[(tID,fID,1)] = 1.0;
            τ2[(fID,tID,1)] = 1.0;
            τ2[(tID,fID,1)] = τ;
            σ[(fID,tID,1)] = σ1/180*pi;
            σ[(tID,fID,1)] = -σ1/180*pi;

          if parse(Float64,brdata[6]) == 0
              rateA[(fID,tID,1)] = Inf;
              rateA[(tID,fID,1)] = Inf;
          else
              rateA[(fID,tID,1)] = parse(Float64,brdata[6])/baseMVA;
              rateA[(tID,fID,1)] = parse(Float64,brdata[6])/baseMVA;
          end

          angmin[(fID,tID,1)] = parse(Float64,brdata[12]);
          angmax[(fID,tID,1)] = parse(Float64,brdata[13]);
          angmin[(tID,fID,1)] = parse(Float64,brdata[12]);
          angmax[(tID,fID,1)] = parse(Float64,brdata[13]);
          push!(brList1,(fID,tID,1));
          push!(brList1,(tID,fID,1));
        else
            m[(fID,tID)] += 1;
            m[(tID,fID)] += 1;
            r = parse(Float64,brdata[3]);
            x = parse(Float64,brdata[4]);
            bc1 = parse(Float64,brdata[5]);
            τ = parse(Float64,brdata[9]);
            if τ == 0.0
                τ = 1.0;
            end
            σ1 = parse(Float64,brdata[10]);
            g[(fID,tID,m[(fID,tID)])] = round(r/(r^2 + x^2),digits=6);
            g[(tID,fID,m[(tID,fID)])] = round(r/(r^2 + x^2),digits=6);
            b[(fID,tID,m[(fID,tID)])] = round(-x/(r^2 + x^2),digits=6);
            b[(tID,fID,m[(tID,fID)])] = round(-x/(r^2 + x^2),digits=6);
            bc[(fID,tID,m[(fID,tID)])] = bc1;
            bc[(tID,fID,m[(tID,fID)])] = bc1;
            τ1[(fID,tID,m[(fID,tID)])] = τ;
            τ1[(tID,fID,m[(tID,fID)])] = 1.0;
            τ2[(fID,tID,m[(tID,fID)])] = 1.0;
            τ2[(tID,fID,m[(tID,fID)])] = τ;
            σ[(fID,tID,m[(fID,tID)])] = σ1/180*pi;
            σ[(tID,fID,m[(tID,fID)])] = -σ1/180*pi;

            if parse(Float64,brdata[6]) == 0
                rateA[(fID,tID,m[(fID,tID)])] = Inf;
                rateA[(tID,fID,m[(tID,fID)])] = Inf;
            else
                rateA[(fID,tID,m[(fID,tID)])] = parse(Float64,brdata[6])/baseMVA;
                rateA[(tID,fID,m[(tID,fID)])] = parse(Float64,brdata[6])/baseMVA;
            end

            angmin[(fID,tID,m[(fID,tID)])] = parse(Float64,brdata[12]);
            angmax[(fID,tID,m[(fID,tID)])] = parse(Float64,brdata[13]);

            angmin[(tID,fID,m[(tID,fID)])] = parse(Float64,brdata[12]);
            angmax[(tID,fID,m[(tID,fID)])] = parse(Float64,brdata[13]);

            push!(brList1,(fID,tID,m[(fID,tID)]));
            push!(brList1,(tID,fID,m[(tID,fID)]));
        end
    end

    # build the reverse dictionary
    brRev = Dict();
    for k in brList1
        brRev[k] = (k[2],k[1],k[3]);
    end

    return brList1,brRev,g,b,bc,angmax,angmin,rateA,τ1,τ2,σ;
end

function parsecST(cST,genIDList,baseMVA)
    # parse the cST string collected from readMP
    # input: the strings containing cost information
    # output: data variable

    # separate the raw strings line by line
    cSTlbl = collect(eachmatch(r"([0-9\t\ \.\-e]+[\;\r]*\n)",cST));
    cp = Dict();
    cq = Dict();

    if length(genIDList) == length(cSTlbl)
        # only active power
        counter = 0;

        for cstr in cSTlbl
            csstr = cstr.captures[1];
            csstr = strip(strip(csstr),';');
            cdata = split(csstr);

            counter += 1;
            typeC = parse(Int64,cdata[1]);
            startC = parse(Float64,cdata[2]);
            endC = parse(Float64,cdata[3]);
            nC = parse(Float64,cdata[4]);
            paramC = [];
            for i in 1:Int64(nC)
                push!(paramC,parse(Float64,cdata[i+4])*(baseMVA^(Int64(nC) - i)));
            end

            citem = costDataType(typeC,startC,endC,nC,paramC);
            cp[genIDList[counter]] = citem;
        end
    elseif length(genIDList)*2 == length(cSTlbl)
        # both active power and reactive power
        counter = 0;
        for cstr in cSTlbl
            csstr = cstr.captures[1];
            csstr = strip(strip(csstr),';');
            cdata = split(csstr,r"\t");

            counter += 1;
            typeC = parse(Int64,cdata[1]);
            startC = parse(Float64,cdata[2]);
            endC = parse(Float64,cdata[3]);
            nC = parse(Float64,cdata[4]);
            paramC = [];
            for i in 1:Int64(nC)
                push!(paramC,parse(Float64,cdata[i+4])*(baseMVA^(Int64(nC) - i)));
            end

            citem = costDataType(typeC,startC,endC,nC,paramC);
            if counter <= length(genIDList)
                cp[genIDList[counter]] = citem;
            else
                cq[genIDList[counter - length(genIDList)]] = citem;
            end
        end
    else
        error("The number of cost entries does not match the number of generators.")
    end

    return cp,cq;
end

function constFixed(baseMVA,bType,IDList,genIDList,brList,brRev,Vmax,Vmin,L,LR,Pmax,Pmin,Qmax,Qmin,gs,bs,Vmag,Vang,Pd,Qd,Pg,Qg,RU,RD,g,b,bc,angmax,angmin,rateA,τ1,τ2,σ,cp,cq,cz)
    # combine all the data to a struct
    connectPair = [];
    connectDict = Dict();
    branchDict1 = Dict();
    branchDict2 = Dict();
    for i in IDList
        connectDict[i] = [];
        branchDict1[i] = [];
        branchDict2[i] = [];
    end
    for k in brList
        push!(branchDict1[k[1]],k);
        push!(branchDict2[k[2]],k);
        if !((k[1],k[2]) in connectPair)
            push!(connectPair,(k[1],k[2]));
            push!(connectDict[k[1]],k[2]);
        end
    end
    kpDict = Dict();
    for k in brList
        if (k[1],k[2]) in keys(kpDict)
            push!(kpDict[(k[1],k[2])],k);
        else
            kpDict[(k[1],k[2])] = [k];
        end
    end

    # collapse all the lines with high conductance
    clusterList = [];
    for k in brList
        if abs(b[k]) >= 1e4
            inCluster = false;
            for item in clusterList
                if (k[1] in item) && !(k[2] in item)
                    push!(item, k[2]);
                    inCluster = true;
                elseif (k[2] in item) && !(k[1] in item)
                    push!(item, k[1]);
                    inCluster = true;
                elseif (k[2] in item) && (k[1] in item)
                    inCluster = true;
                end
            end
            if !inCluster
                push!(clusterList,[k[1],k[2]]);
            end
        end
    end

    # for each cluster, combine the nodes to the lowest indexed node
    for item in clusterList
        # combine the nodes to the lowest indexed node
        combinedN = minimum(item);
        item1 = copy(item);
        filter!(e -> e != combinedN, item1);
        # remove the node
        filter!(e -> !(e in item1), IDList);
        for i in item1
            # remove the pseudo line
            filter!(e -> !(((e[1] == i)&&(e[2] == combinedN))||((e[1] == combinedN)&&(e[2] == i))), brList);
            for k in keys(brRev)
                if !(k in brList)
                    delete!(brRev,k);
                    delete!(kpDict,(k[1],k[2]));
                    filter!(e -> e!=k, branchDict1[i]);
                    filter!(e -> e!=k, branchDict2[i]);
                    delete!(rateA,k);
                    delete!(g,k);
                    delete!(b,k);
                    delete!(bc,k);
                    delete!(τ1,k);
                    delete!(τ2,k);
                    delete!(σ,k);
                end
            end
            # consolidate the generator information
            if i in keys(LR)
                rmGenList = LR[i];
                if combinedN in keys(LR)
                    append!(LR[combinedN], rmGenList);
                    for gitem in rmGenList
                        L[gitem] = combinedN;
                    end
                else
                    LR[combinedN] = rmGenList;
                    for gitem in rmGenList
                        L[gitem] = combinedN;
                    end
                end
                delete!(LR,i);
            end

            # remove the node information
            Vmax[combinedN] = min(Vmax[combinedN],Vmax[i]);
            Vmin[combinedN] = max(Vmin[i],Vmin[combinedN]);
            gs[combinedN] += gs[i];
            bs[combinedN] += bs[i];
            Pd[combinedN] += Pd[i];
            Qd[combinedN] += Qd[i];
            delete!(Vmax,i);
            delete!(Vmin,i);
            delete!(gs,i);
            delete!(bs,i);
            delete!(Pd,i);
            delete!(Qd,i);
            delete!(Vmag,i);
            delete!(Vang,i);
            delete!(angmax,i);
            delete!(angmin,i);

            # consolidate the arcs
            branchDict1Temp = [];
            for k in branchDict1[i]
                filter!(e -> e != k, brList);
                if (combinedN,k[2]) in keys(kpDict)
                    newbr = (combinedN,k[2],length(kpDict[(combinedN,k[2])]) + 1);
                    push!(brList,newbr);
                    brRev[newbr] = (newbr[2],newbr[1],newbr[3]);
                    kpDict[(newbr[1],newbr[2])] = newbr;
                else
                    newbr = (combinedN,k[2],1);
                    push!(brList,newbr);
                    brRev[newbr] = (newbr[2],newbr[1],newbr[3]);
                    kpDict[(newbr[1],newbr[2])] = newbr;
                end
                filter!(e -> e != k, branchDict2[k[2]]);
                push!(branchDict2[k[2]],newbr);
                rateA[newbr] = rateA[k];
                g[newbr] = g[k];
                b[newbr] = b[k];
                bc[newbr] = bc[k];
                τ1[newbr] = τ1[k];
                τ2[newbr] = τ2[k];
                σ[newbr] = σ[k];
                push!(branchDict1Temp,newbr);
                delete!(brRev,k);
                delete!(kpDict,(k[1],k[2]));
                delete!(rateA,k);
                delete!(g,k);
                delete!(b,k);
                delete!(bc,k);
                delete!(τ1,k);
                delete!(τ2,k);
                delete!(σ,k);
            end
            branchDict1[i] = branchDict1Temp;
            branchDict2Temp = [];
            for k in branchDict2[i]
                filter!(e -> e != k, brList);
                if (k[1],combinedN) in keys(kpDict)
                    newbr = (k[1],combinedN,length(kpDict[(k[1],combinedN)]) + 1);
                    push!(brList,newbr);
                    brRev[newbr] = (newbr[2],newbr[1],newbr[3]);
                    kpDict[(newbr[1],newbr[2])] = newbr;
                else
                    newbr = (k[1],combinedN,1);
                    push!(brList,newbr);
                    brRev[newbr] = (newbr[2],newbr[1],newbr[3]);
                    kpDict[(newbr[1],newbr[2])] = newbr;
                end
                filter!(e -> e != k, branchDict1[k[1]]);
                push!(branchDict1[k[1]],newbr);
                rateA[newbr] = rateA[k];
                g[newbr] = g[k];
                b[newbr] = b[k];
                bc[newbr] = bc[k];
                τ1[newbr] = τ1[k];
                τ2[newbr] = τ2[k];
                σ[newbr] = σ[k];
                push!(branchDict2Temp,newbr);
                delete!(brRev,k);
                delete!(kpDict,(k[1],k[2]));
                delete!(rateA,k);
                delete!(g,k);
                delete!(b,k);
                delete!(bc,k);
                delete!(τ1,k);
                delete!(τ2,k);
                delete!(σ,k);
            end
            branchDict2[i] = branchDict2Temp;
        end
    end

    ii = 0;
    busInd = Dict();
    for i in IDList
        ii += 1;
        busInd[i] = ii;
    end

    # obtain the pairs that are connected
    connectPair = [];
    connectDict = Dict();
    branchDict1 = Dict();
    branchDict2 = Dict();
    for i in IDList
        connectDict[i] = [];
        branchDict1[i] = [];
        branchDict2[i] = [];
    end
    for k in brList
        push!(branchDict1[k[1]],k);
        push!(branchDict2[k[2]],k);
        if !((k[1],k[2]) in connectPair)
            push!(connectPair,(k[1],k[2]));
            push!(connectDict[k[1]],k[2]);
        end
    end
    kpDict = Dict();
    for k in brList
        if (k[1],k[2]) in keys(kpDict)
            push!(kpDict[(k[1],k[2])],k);
        else
            kpDict[(k[1],k[2])] = [k];
        end
    end
    RUadj = Dict();
    RDadj = Dict();
    for iKey in keys(RU)
        RUadj[iKey] = RU[iKey]/baseMVA;
    end
    for iKey in keys(RD)
        RDadj[iKey] = RD[iKey]/baseMVA;
    end
    fData = fixedData(baseMVA,bType,IDList,genIDList,brList,brRev,
                        L,LR,Vmax,Vmin,Pmax,Pmin,Qmax,Qmin,gs,bs,Vmag,Vang,
                        Pd,Qd,Pg,Qg,RUadj,RDadj,g,b,bc,angmax,angmin,rateA,τ1,τ2,σ,cp,cq,cz,
                        busInd,branchDict1,branchDict2,connectPair,connectDict,kpDict);
    return fData
end

function readStatic(fileName::String, cz = 1e7)
    # read in the static network information from a .m file (MatPower)
    busST,genST,brST,cST,uST,baseMVA = readMP(fileName);
    IDList,Vmax,Vmin,gs,bs,Vmag,Vang,Pd,Qd,bType = parsebusST(busST,baseMVA);
    genIDList,Loc,LocRev,Pmax,Pmin,Qmax,Qmin,Pg,Qg,RU,RD = parsegenST(genST,baseMVA);
    brList,brRev,g,b,bc,angmax,angmin,rateA,τ1,τ2,σ = parsebrST(brST,baseMVA);
    cp,cq = parsecST(cST,genIDList,baseMVA);
    fData = constFixed(baseMVA,bType,IDList,genIDList,brList,brRev,Vmax,Vmin,Loc,LocRev,Pmax,Pmin,Qmax,Qmin,
            gs,bs,Vmag,Vang,Pd,Qd,Pg,Qg,RU,RD,g,b,bc,angmax,angmin,rateA,τ1,τ2,σ,cp,cq,cz);

    # specify the line power flow constraint if there is none
    θu = pi/3;
    for k in fData.brList
        if fData.rateA[k] == Inf
            fData.rateA[k] = sqrt(fData.g[k]^2+fData.b[k]^2)*max(fData.Vmax[k[1]],fData.Vmax[k[2]])*sqrt(fData.Vmax[k[1]]^2 + fData.Vmax[k[2]]^2 - 2*fData.Vmax[k[1]]*fData.Vmin[k[2]]*cos(θu));
        end
    end
    return fData;
end

function readDisruption(tfileName,ofileName,fileType,τ = nothing)
    # read in the disruption data: disruption time distribution
    if fileType == "csv"
        # csv file format:
        # First row: time period t
        # Second row: pr{D = t}
        # Third row: Ω
        # Fourth row: pr{B_ω = 0}
        dataT = readdlm(tfileName, ',');
        tDistrn = Dict();
        mt,nt = size(dataT);
        for tInd in 1:nt
            tDistrn[Int64(dataT[1,tInd])] = dataT[2,tInd];
        end

        # detect the dimension of Ω
        dataOme = readdlm(ofileName, ',');
        ωDistrn = Dict();
        ωDict = Dict();
        ωτ = Dict();
        mo,no = size(dataOme);
        for oInd in 1:mo
            scenInd = dataOme[oInd,1];
            ωDistrn[scenInd] = dataOme[oInd,2];
            ωDict[scenInd] = [];
            if τ == nothing
                τr = dataOme[oInd,3];
            else
                τr = τ;
            end
            ωτ[scenInd] = τr;
            for noInd in 4:no
                lineStr = dataOme[oInd,noInd];
                if ',' in lineStr
                    # parse the line
                    startN,endN = match(r"\(([0-9]+),([0-9]+)\)",lineStr).captures;
                    push!(ωDict[scenInd],(parse(Int64,startN),parse(Int64,endN)));
                elseif !(lineStr == "")
                    # parse the node
                    disN = lineStr;
                    push!(ωDict[scenInd],disN);
                end
            end
        end
        if sum(values(ωDistrn)) == 1
            pDistr = probDistrn(tDistrn,ωDistrn,ωDict,ωτ);
            return pDistr;
        else
            ωList = [item for item in keys(ωDistrn)];
            ωDistrn[ωList[1]] = 1 - sum([ωDistrn[item] for item in keys(ωDistrn) if item != ωList[1]]);
            pDistr = probDistrn(tDistrn,ωDistrn,ωDict,ωτ);
            return pDistr;
        end
    else
        println("Currently your file type is not supported");
    end
end

function readDisruption_old(fileName,fileType)
    # read in the disruption data: disruption time distribution
    if fileType == "csv"
        # csv file format:
        # First row: time period t
        # Second row: pr{D = t}
        # Third row: Ω
        # Fourth row: pr{B_ω = 0}
        dataRaw = readdlm(fileName, ',');

        # detect the dimension of T
        dataT = dataRaw[1:2,:];
        emptyInd = findfirst(x -> x=="", dataT[1,:]);
        if emptyInd != nothing
            # if it is not full length
            dataT = dataT[:,1:(emptyInd - 1)];
        end
        tDistrn = Dict();
        mt,nt = size(dataT);
        for tInd in 1:nt
            tDistrn[Int64(dataT[1,tInd])] = dataT[2,tInd];
        end

        # detect the dimension of Ω
        dataOme = dataRaw[3:4,:];
        emptyInd = findfirst(x -> x=="", dataOme[1,:]);
        if emptyInd != nothing
            # if it is not full length
            dataOme = dataOme[:,1:(emptyInd - 1)];
        end
        ωDistrn = Dict();
        mo,no = size(dataOme);
        for oInd in 1:no
            lineStr = dataOme[1,oInd];
            if ',' in lineStr
                # parse the line
                startN,endN = match(r"\(([0-9]+),([0-9]+)\)",lineStr).captures;
                ωDistrn[parse(Int64,startN),parse(Int64,endN)] = dataOme[2,oInd];
            else
                # parse the node
                disN = lineStr;
                ωDistrn[disN] = dataOme[2,oInd];
            end
        end
        if sum(values(ωDistrn)) == 1
            pDistr = probDistrn(tDistrn,ωDistrn);
            return pDistr;
        else
            ωList = [item for item in keys(ωDistrn)];
            ωDistrn[ωList[1]] = 1 - sum([ωDistrn[item] for item in keys(ωDistrn) if item != ωList[1]]);
            pDistr = probDistrn(tDistrn,ωDistrn);
            return pDistr;
        end
    else
        println("Currently your file type is not supported");
    end
end

function readBattery(fileName,fileType,baseMVA = 100)
    # read in the battery information: charging/discharging factor, capacity
    if fileType == "csv"
        # csv file format:
        # First column: node ID
        # Second column: ηd
        # Third column: ηc
        # Fourth column: capacity
        # Fifth column: cost (optional)
        dataRaw = readdlm(fileName, ',');
        mb,nb = size(dataRaw);
        capacity = Dict();
        cost = Dict();
        ηα = Dict();
        ηβ = Dict();
        IDList = [];
        LocDict = Dict();
        bInv = Dict();
        uCap = Dict();

        for i in 1:mb
            ID = Int64(dataRaw[i,1]);
            push!(IDList,ID);
            loc = Int64(dataRaw[i,2]);
            LocDict[ID] = loc;
            capacity[ID] = dataRaw[i,3]/baseMVA;
            cost[ID] = dataRaw[i,4]*baseMVA;
            bInv[ID] = dataRaw[i,5]/baseMVA;
            uCap[ID] = dataRaw[i,6]/baseMVA;
            ηparams = [j for j in dataRaw[i,7:nb] if j != ""];
            ηα[ID] = [];
            ηβ[ID] = [];
            for j in 1:2:length(ηparams)
                push!(ηα[ID],ηparams[j]);
                # ηα must be monotonously decreasing
                push!(ηβ[ID],ηparams[j + 1]);
            end
        end
        bData = batteryData(IDList,LocDict,bInv,ηα,ηβ,capacity,cost,uCap);
    else
        println("Currently your file type is not supported");
    end
    return bData;
end

function readDemand(fileNameP, fileNameQ, fileType)
    # read in the battery information: charging/discharging factor, capacity
    if fileType == "csv"
        # csv file format:
        # First column: node ID
        # Second column: ηd
        # Third column: ηc
        dataRawP = readdlm(fileNameP, ',');
        mp,np = size(dataRawP);
        pd = Dict();
        for i in 1:mp
            pd[Int64(dataRawP[i,1])] = dataRawP[i,2:np];
        end

        dataRawQ = readdlm(fileNameQ, ',');
        mq,nq = size(dataRawQ);
        qd = Dict();
        for i in 1:mq
            qd[Int64(dataRawQ[i,1])] = dataRawQ[i,2:nq];
        end
        T = nq - 1;
        dData = demandData(T,pd,qd);
        return dData;
    else
        println("Currently your file type is not supported");
    end
end

function readInData(i,caseList,T,τ,cz = 1e4,λD = 0)
    fileAdd = "case$(caseList[i])_ieee.m";
    global fData = readStatic(fileAdd,cz);
    disAddt = "testProbReadt_$(caseList[i]).csv";
    disAddo = "testProbReado_$(caseList[i]).csv";
    global pDistr = readDisruption(disAddt,disAddo,"csv",τ);
    pAdd = "testDataP_$(caseList[i]).csv";
    qAdd = "testDataQ_$(caseList[i]).csv";
    global dData = readDemand(pAdd,qAdd,"csv");
    bAdd = "testDataB_$(caseList[i]).csv";
    global bData = readBattery(bAdd,"csv",fData.baseMVA);

    if λD == 0
        global pDistr = modifyT(pDistr,4/T,T);
    else
        global pDistr = modifyT(pDistr,λD,T);
    end
end

function readInData_tau(i,caseList,T,cz = 1e4,λD = 0,suffix = "",τ = nothing)
    fileAdd = "case$(caseList[i])_ieee.m";
    global fData = readStatic(fileAdd,cz);
    disAddt = "testProbReadt_$(caseList[i]).csv";
    disAddo = "testProbReado_$(caseList[i])_tau$(suffix).csv";
    global pDistr = readDisruption(disAddt,disAddo,"csv",τ);
    pAdd = "testDataP_$(caseList[i]).csv";
    qAdd = "testDataQ_$(caseList[i]).csv";
    global dData = readDemand(pAdd,qAdd,"csv");
    bAdd = "testDataB_$(caseList[i]).csv";
    global bData = readBattery(bAdd,"csv",fData.baseMVA);

    if λD == 0
        global pDistr = modifyT(pDistr,4/T,T);
    else
        global pDistr = modifyT(pDistr,λD,T);
    end
end

function readInData_old(T,ωSet,cz = 500,λD = 1/2)
    fileAdd = "case13_ieee_old.m";
    global fData = readStatic(fileAdd,cz);
    disAdd = "testProbRead_96.csv"
    global pDistr = readDisruption(disAdd,"csv");
    pAdd = "testDataP_96_Fixed.csv";
    qAdd = "testDataQ_96_Fixed.csv";
    global dData = readDemand(pAdd,qAdd,"csv");
    bAdd = "testDataB.csv";
    global bData = readBattery(bAdd,"csv",fData.baseMVA);

    for ω in keys(pDistr.ωDistrn)
        if !(ω in ωSet)
            global pDistr = modifyOmega(pDistr,ω);
        end
    end

    if λD == 0
        global pDistr = modifyT(pDistr,4/T,T);
    else
        global pDistr = modifyT(pDistr,λD,T);
    end
end
