function project570()
    clear all;
    tic
    
    % NODE Id's
    A=1;
    B=2;
    
    C=3;
    D=4;
    
    % Matrix init for routing
    M = zeros(4,4);
    % Reference:
    % http://en.wikipedia.org/wiki/Discrete_event_simulation#Events_List
    % Events List: The simulation maintains at least one list of simulation
    % events. This is sometimes called the pending event set because it
    % lists events that are pending as a result of previously simulated
    % event but have yet to be simulated themselves.
    % I have a elist as an (n by 6) matrix where n is the number of packets
    % waiting in the queue of nodes A and B. We define a packet to be a
    % row in the elist and is of the form
    % [SRC GENTIME TXTIME RXTIME CURTIME COLLISIONS]
    elist=[];
    elist1=[];
    
    % Columns in the event list
    SRC=1;      % Id of the source node
    GENTIME=2;  % Time at which the packet is created at the source node
    TXTIME=3;   % Time at which the packet is transmitted
    RXTIME=4;   % Time at which the packet is received at the dest node
    % The current simulation time w.r.t the packet. You can think of
    % CURTIME to be the time at which the simulator had last visited this
    % packet.
    CURTIME=5;  
    COLLISIONS=6; % Number of collisions incurred by 
                    % this packet due to channel contention
    DESTINATION=7;                
    ROUTINGPATH=8;
    PROPDELAY=9;
    % Reference:
    % http://en.wikipedia.org/wiki/Discrete_event_simulation#Clock 
    % Clock:
    % The simulation must keep track of the current simulation time,
    % whatever measurement units are suitable for the system being modeled.
    % In discrete-event simulations, as opposed to real time simulations,
    % time ‘hops’ because events are instantaneous – the clock skips to the
    % next event start time as the simulation proceeds. 
    CLOCK=0;
    CLOCK1=0;
    
    % Initialization
    frameSize = 1000*8;
    LANTransRate = 100;
    distance = 2000;
    propSpeed = 200;  % 2 x 10^8
    
    TOTALSIM=60*10^4; % Total simulation time
    lambda = .5;
    frameslot = 50; % frame slot time (usec) <-- changed from 500us
    td = frameSize/LANTransRate; %80;% transmission delay on BUS (usec)
    pd = distance/propSpeed; %10; % propagation delay on BUS
    tdelay=td+pd; % total delay incurred during a pkt transmission
	tbackoff = frameslot; % time slot (usec) for backoff algorithm
    maxbackoff = 3; % maximum backoff time is 2^3 frame slot
	
    RLTransRate = 1000;
    % Routing params
    tdr = frameSize/RLTransRate; %8; Transmission delay in the router links
    pdr = distance/propSpeed; %10; Propagation delay in the router links
    tdelayr = tdr + pdr;
    
    LinkUpdateInterval = 2000; % Time interval between change in router link cost
    % the time at which the last packet was generated at node A and B resp.
    GENTIMECURSOR=[0 0 0 0];
        
    % Create packet for node A. Function createpacket() is defined below.
    % check the function first for better understanding.
    createpacket(A); 
    createpacket(B);
    createpacket1(C);
    createpacket1(D);
  
    % Set the CLOCK time. Look at the function for more details. 
    updateclock();
    
    if size(elist,1) == 0
        disp('No packets to simulate');
        return;
    end
    
    if size(elist1,1) == 0
        disp('No packets to simulate');
        return;
    end
    
    
    % Collect the statistics in this array.
    SIMRESULT=[];
    SIMRESULT1=[];
        
    while(1)
        % Update the clock.
        updateclock();     
        % Find the source node of the packet at the first row in elist.
        src=elist(1,SRC);
        secsrc=1;
        if(src == 1)
            secsrc=2;
        end 
        
        src1=elist1(1,SRC);
        secsrc1 = 3;
        if(src1 == 3)
            secsrc1 = 4;
        end    
        
		%dst=mod(src,2)+1;
        timediff=elist(2,CURTIME)-elist(1,CURTIME);
        timediff1=elist1(2,CURTIME)-elist1(1,CURTIME);
        
        dst = elist(1, DESTINATION);
        dst1 = elist1(1, DESTINATION);
        
		if timediff > pd
            % routing
            SRSize = size(SIMRESULT);
            PrevTxTime = 0;
            SRLastRow = SRSize(1,1);
            if SRLastRow > 0
                PrevTxTime = SIMRESULT(SRLastRow, 3);
            end
            updateRouteFlag = 1;
            timegap = elist(1,CURTIME) - PrevTxTime;
            if timegap < LinkUpdateInterval && PrevTxTime > 0
                updateRouteFlag = 0; 
            end
            if dst == 3  || dst == 4
                routingdelay = getRoutingPath(src, updateRouteFlag);
            end
			%% No collision case
            % Set the tx time. This should be the time when the packet
            % is transmitted for the first time.
            if elist(1,TXTIME)==0
                elist(1,TXTIME)=elist(1,CURTIME);
            end
            % Set the rx time.
            elist(1,RXTIME)=elist(1,CURTIME)+tdelay;
            if dst == 3  || dst == 4
                elist(1,RXTIME)=elist(1,CURTIME)+(2 * tdelay) + routingdelay;
            else
                elist(1, PROPDELAY) = tdelay; 
            end
            updatesimlist();
			
            createpacket(src);
			
			% add tdelay to CLOCK. Check delaypkts() for more details.
            delaypkts(tdelay, 1);
        else
            % Collision and so Backoff
            % Collision is detected by src node after 2*(timediff
            % +(propdelay-timediff)/2)=timediff+propdelay
            % Collision is detected by interfering node after 
            % 2*(propdelay-timediff)/2=propdelay-timediff
            %disp('backoff');
            if elist(1,TXTIME)==0
                elist(1,TXTIME)=elist(1,CURTIME);
            end
            if elist(2,TXTIME)==0
                elist(2,TXTIME)=elist(2,CURTIME);
            end
            elist(1,COLLISIONS)=elist(1,COLLISIONS)+1;
            elist(2,COLLISIONS)=elist(2,COLLISIONS)+1;
            
            if (elist(1,COLLISIONS)<maxbackoff)
                bk(src)=(randi(2^(elist(1,COLLISIONS)),1,1)-1)*tbackoff;
            else
                bk(src)=(randi(2^(maxbackoff),1,1)-1)*tbackoff;
            end
            if (elist(2,COLLISIONS)<maxbackoff)
                bk(secsrc)=(randi(2^(elist(2,COLLISIONS)),1,1)-1)*tbackoff;
            else
                bk(secsrc)=(randi(2^(maxbackoff),1,1)-1)*tbackoff;
            end
			
			% Delay the packets at each node.
            delaynodepkts(src,pd+timediff+bk(src),1);
            delaynodepkts(secsrc,pd-timediff+bk(secsrc),1);
        end
		
        
        %********************
        
        if timediff1 > pd
            % routing
            SRSize1 = size(SIMRESULT1);
            PrevTxTime1 = 0;
            SRLastRow1 = SRSize1(1,1);
            if SRLastRow1 > 0
                PrevTxTime1 = SIMRESULT1(SRLastRow1, 3);
            end
            updateRouteFlag1 = 1;
            timegap1 = elist1(1,CURTIME) - PrevTxTime1;
            if timegap1 < LinkUpdateInterval && PrevTxTime1 > 0
                updateRouteFlag1 = 0;
            end
            if dst1 == 1  || dst1 == 2
                routingdelay = getRoutingPath1(src1, updateRouteFlag1);
            else
                elist1(1, PROPDELAY) = tdelay; 
            end
			%% No collision case
            % Set the tx time. This should be the time when the packet
            % is transmitted for the first time.
            if elist1(1,TXTIME)==0
                elist1(1,TXTIME)=elist1(1,CURTIME);
            end
            % Set the rx time.
            elist1(1,RXTIME)=elist1(1,CURTIME)+tdelay;
            if dst1 == 1  || dst1 == 2
                elist1(1,RXTIME)=elist1(1,CURTIME)+ (2 * tdelay) + routingdelay;
            end
            updatesimlist1();
			
            createpacket1(src1);
			
			% add tdelay to CLOCK. Check delaypkts() for more details.
            delaypkts(tdelay, 2);
        else
            % Collision and so Backoff
            % Collision is detected by src node after 2*(timediff
            % +(propdelay-timediff)/2)=timediff+propdelay
            % Collision is detected by interfering node after 
            % 2*(propdelay-timediff)/2=propdelay-timediff
            %disp('backoff');
            if elist1(1,TXTIME)==0
                elist1(1,TXTIME)=elist1(1,CURTIME);
            end
            if elist1(2,TXTIME)==0
                elist1(2,TXTIME)=elist1(2,CURTIME);
            end
            elist1(1,COLLISIONS)=elist1(1,COLLISIONS)+1;
            elist1(2,COLLISIONS)=elist1(2,COLLISIONS)+1;
            
            if (elist1(1,COLLISIONS)<maxbackoff)
                bk(src1)=(randi(2^(elist1(1,COLLISIONS)),1,1)-1)*tbackoff;
            else
                bk(src1)=(randi(2^(maxbackoff),1,1)-1)*tbackoff;
            end
            if (elist1(2,COLLISIONS)<maxbackoff)
                bk(secsrc1)=(randi(2^(elist1(2,COLLISIONS)),1,1)-1)*tbackoff;
            else
                bk(secsrc1)=(randi(2^(maxbackoff),1,1)-1)*tbackoff;
            end
			
			% Delay the packets at each node.
            delaynodepkts(src1,pd+timediff+bk(src1),2);
            delaynodepkts(secsrc1,pd-timediff+bk(secsrc1),2);
        end
        
        %********************
        
        
        if min(GENTIMECURSOR) > TOTALSIM
            disp('Completed!');
            disp(SIMRESULT);
            disp(SIMRESULT1);
            calcstat();
            break;
        end
    end
       
	   function calcstat()
        FromA=SIMRESULT(SIMRESULT(:,SRC)==A,:);
        FromB=SIMRESULT(SIMRESULT(:,SRC)==B,:);
        FromC=SIMRESULT1(SIMRESULT1(:,SRC)==C,:);
        FromD=SIMRESULT1(SIMRESULT1(:,SRC)==D,:);
        
        % Total packets sent
        Anum=length(FromA);
        Bnum=length(FromB);
        Cnum=length(FromC);
        Dnum=length(FromD);
        
        % Queue delay
        queuedelaya=FromA(:,RXTIME)-FromA(:,GENTIME);
        queuedelayb=FromB(:,RXTIME)-FromB(:,GENTIME);
        queuedelayc=FromC(:,RXTIME)-FromC(:,GENTIME);
        queuedelayd=FromD(:,RXTIME)-FromD(:,GENTIME);
        
        queuedelaya=queuedelaya-tdelay;
        queuedelayb=queuedelayb-tdelay;
        queuedelayc=queuedelayc-tdelay;
        queuedelayd=queuedelayd-tdelay;
        
        
        % Access delay
        accessdelaya=FromA(:,RXTIME)-FromA(:,TXTIME);
        accessdelayb=FromB(:,RXTIME)-FromB(:,TXTIME);
        
        figure;
        subplot(3,1,1)
        plot(1:size(queuedelaya,1),queuedelaya(1:end))
        axis([0 Anum 0 max(queuedelaya)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Queue delay at node A');
         
              
        subplot(3,1,2)
        plot(1:size(accessdelaya,1),accessdelaya(1:end))
        axis([0 Anum 0 max(accessdelaya)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Access delay at node A');
        
        
         % Frame interval
        subplot(3,1,3)
        frameinta=FromA(2:end,GENTIME)-FromA(1:end-1,GENTIME);
        plot(1:Anum-1,frameinta/1000)
        axis([0 Anum 0 max(frameinta)/1000])
        xlabel('Packet sequence #');
        ylabel('Frame interval in msec');
        title('Frame intervals at node A');
        
        
        figure;
        subplot(3,1,1)
        plot(1:size(queuedelayb,1),queuedelayb(1:end))
        axis([0 Bnum 0 max(queuedelayb)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Queue delay at node B');
           
              
        subplot(3,1,2)
        plot(1:size(accessdelayb,1),accessdelayb(1:end))
        axis([0 Bnum 0 max(accessdelayb)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Access delay at node B');
        
       
        
        subplot(3,1,3)
        frameintb=FromB(2:end,GENTIME)-FromB(1:end-1,GENTIME);
        plot(1:Bnum-1,frameintb/1000)
        axis([0 Bnum 0 max(frameintb)/1000])
        xlabel('Packet sequence #');
        ylabel('Frame interval in msec');
        title('Frame intervals at node B');
        
        % Access delay
        accessdelayc=FromC(:,RXTIME)-FromC(:,TXTIME);
        accessdelayd=FromD(:,RXTIME)-FromD(:,TXTIME);
        
        figure;
        subplot(3,1,1)
        plot(1:size(queuedelayc,1),queuedelayc(1:end))
        axis([0 Cnum 0 max(queuedelayc)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Queue delay at node C');
        
        
        subplot(3,1,2)
        plot(1:size(accessdelayc,1),accessdelayc(1:end))
        axis([0 Cnum 0 max(accessdelayc)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Access delay at node C');
                   
        
        % Frame interval
        subplot(3,1,3)
        frameintc=FromC(2:end,GENTIME)-FromC(1:end-1,GENTIME);
        plot(1:Cnum-1,frameintc/1000)
        axis([0 Cnum 0 max(frameintc)/1000])
        xlabel('Packet sequence #');
        ylabel('Frame interval in msec');
        title('Frame intervals at node C');
        
        
        
        figure;
        
        subplot(3,1,1)
        plot(1:size(queuedelayd,1),queuedelayd(1:end))
        axis([0 Dnum 0 max(queuedelayd)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Queue delay at node D');
        
        
        subplot(3,1,2)
        plot(1:size(accessdelayd,1),accessdelayd(1:end))
        axis([0 Dnum 0 max(accessdelayd)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Access delay at node D');
        
        
        
        subplot(3,1,3)
        frameintd=FromD(2:end,GENTIME)-FromD(1:end-1,GENTIME);
        plot(1:Dnum-1,frameintd/1000)
        axis([0 Dnum 0 max(frameintd)/1000])
        xlabel('Packet sequence #');
        ylabel('Frame interval in msec');
        title('Frame intervals at node D');
        
        
        
        
        
        
        
        % Histogram to verify if the distribution is exp
        figure;
        subplot(2,2,1)
        hist(frameinta,60);
        xlabel('frame intervals in \mu sec');
        ylabel('# of frames');
        subplot(2,2,2)
        hist(frameintb,60);
        xlabel('frame intervals in \mu sec');
        ylabel('# of frames');
        subplot(2,2,3)
        hist(frameintc,60);
        xlabel('frame intervals in \mu sec');
        ylabel('# of frames');
        subplot(2,2,4)
        hist(frameintd,60);
        xlabel('frame intervals in \mu sec');
        ylabel('# of frames');
        
        % Calculate count of paths
        plotRoutingCounts();
        
        % Calculate propagation delay
        plotPropDelayFromAandB();
        plotPropDelayFromCandD();
        
        % mean end to end time
        meanendtoenda=mean(FromA(:,RXTIME)-FromA(:,GENTIME));
        
        meanendtoendb=mean(FromB(:,RXTIME)-FromB(:,GENTIME));
        meanendtoendc=mean(FromC(:,RXTIME)-FromC(:,GENTIME));
        
        meanendtoendd=mean(FromD(:,RXTIME)-FromD(:,GENTIME));
        
        % Average end to end throughput
        Avgtha=((1000*8)/meanendtoenda)*10^6;  % bits/sec
        Avgthb=((1000*8)/meanendtoendb)*10^6;  % bits/sec
        Avgthc=((1000*8)/meanendtoendc)*10^6;  % bits/sec
        Avgthd=((1000*8)/meanendtoendd)*10^6;  % bits/sec
        avgth2=mean(8000.*(FromA(:,RXTIME)-FromA(:,GENTIME)).^(-1));
        
        fprintf('Total packets sent from node A=%d\n',Anum);
        fprintf('Total packets sent from node B=%d\n',Bnum);
        fprintf('Total packets sent from node C=%d\n',Cnum);
        fprintf('Total packets sent from node D=%d\n',Dnum);
        fprintf('Average frame interval at node A=%d\n',mean(frameinta));
        fprintf('Average frame interval at node B=%d\n',mean(frameintb));
        fprintf('Average frame interval at node C=%d\n',mean(frameintc));
        fprintf('Average frame interval at node D=%d\n',mean(frameintd));
        fprintf('Average access delay at node A=%d\n',mean(accessdelaya));
        fprintf('Average access delay at node B=%d\n',mean(accessdelayb));
        fprintf('Average access delay at node C=%d\n',mean(accessdelayc));
        fprintf('Average access delay at node D=%d\n',mean(accessdelayd));
        fprintf('Average queue delay at node A=%d\n',mean(queuedelaya));
        fprintf('Average queue delay at node B=%d\n',mean(queuedelayb));
        fprintf('Average queue delay at node C=%d\n',mean(queuedelayc));
        fprintf('Average queue delay at node D=%d\n',mean(queuedelayd));
        fprintf('Average end to end throughput from A to B=%d\n',Avgtha);
        fprintf('Average end to end throughput from B to A=%d\n',Avgthb);
        fprintf('Average end to end throughput from A to B=%d\n',Avgthc);
        fprintf('Average end to end throughput from B to A=%d\n',Avgthd);
        
        fprintf('Simulation end time=%d\n',CLOCK);
       % fprintf('Simulation end time=%d\n',CLOCK1);
        
    end
	
    % The clock slips to the RXTIME i.e., add delay time to CLOCK.
    function delaypkts(delay,pkt)
        if(pkt == 1)
            CLOCK=CLOCK+delay;
            list=find(((elist(:,CURTIME)-CLOCK) < 0));
            elist(list,CURTIME)=CLOCK;
        else        
            CLOCK1=CLOCK1+delay;
            list=find(((elist1(:,CURTIME)-CLOCK) < 0));
            elist1(list,CURTIME)=CLOCK;
        end
        % It might so happen that the new packet at the SRC node might have
        % been generated when the previous packet was in flight. This new
        % packet cannot be transmitted immediately and hence has to wait
        % till the previos packet has reached the destination. 
        
        % list will have the row number of elist whose CURTIME field value
        % is less than CLOCK. Remember CLOCK is now the RXTIME. 
            %list=find(((elist(:,CURTIME)-CLOCK) < 0));
        % Set the CURTIME field of all the rows in list to CLOCK.
            %elist(list,CURTIME)=CLOCK;
    end
   
    function updateclock()
        % SORTROWS(elist,CURTIME) sorts the rows of elist in ascending
        % order for the column CURTIME.
        elist=sortrows(elist,CURTIME);
        elist1=sortrows(elist1,CURTIME);
        % Set the clock to the CURTIME of the packet in the first row of
        % the elist since this is the packet that contends first for the
        % channel. 
        CLOCK=elist(1,CURTIME);
        CLOCK1=elist1(1,CURTIME);
    end
    
    % If the number of arrivals in any given time interval [0,t] follows
    % the Poisson distribution, with mean = \lambda \cdot t, then the lengths of the
    % inter-arrival times follow the Exponential distribution, with mean
    % 1/\lambda. We can create packets all at once, i.e, keep creating packets
    % until the cummulative sum of the inter-arrival time is greater than
    % TOTALSIM (60 seconds). You might have around 10^4 packets generated
    % and all these packets will go into the event list. However if think
    % carefully, only two packets contend for the channel access at any
    % time, these are the packets at the head of the queue in node A and
    % node B. So you can avoid having other packets in the event list. 
    % Hence a smarter way is to create a new packet at a node soon
    % after the packet at the head of the queue has left the node. However
    % you must be care keep track of the cummulative sum of inter-arrival
    % time upto the previous packet (the packet that just left the queue)
    % inorder to calculate the birth time of the packet generated. 
    
    % Note: You can either create all packets at once or follow the
    % approach I have taken. It is totally upto you. The outcome will be
    % same with both the approaches will lead to the same outcome.
    function pkt=createpacket(nodeid)
        % Find the inter-arrival time.
        interarvtime = round(frameslot*exprnd(1/lambda,1,1));
        
        % Find the birth time.
        GENTIMECURSOR(nodeid)=GENTIMECURSOR(nodeid)+interarvtime;
        
        % Create the packet. Unknown fields are set to 0.
        % [SRC= nodeid GENTIME=birthtime TXTIME=0, RXTIME=0
        % CURTIME=birthtime COLLISIONS=0]
        curdest = randi(4,1,1);
        while(curdest == nodeid)
            curdest = randi(4,1,1);
        end
        
        pkt=[nodeid GENTIMECURSOR(nodeid) 0 0 GENTIMECURSOR(nodeid) 0 curdest 0 0];
        
        % Enqueue to the event list. In matlab you can append a row x to a
        % matrix X by using the command X=[X; x]
        elist=[elist; pkt];
    end
	
    function pkt=createpacket1(nodeid)
        % Find the inter-arrival time.
        interarvtime = round(frameslot*exprnd(1/lambda,1,1));
        
        % Find the birth time.
        GENTIMECURSOR(nodeid)=GENTIMECURSOR(nodeid)+interarvtime;
        
        % Create the packet. Unknown fields are set to 0.
        % [SRC= nodeid GENTIME=birthtime TXTIME=0, RXTIME=0
        % CURTIME=birthtime COLLISIONS=0]
        curdest = randi(4,1,1);
        while(curdest == nodeid)
            curdest = randi(4,1,1);
        end
        
        pkt=[nodeid GENTIMECURSOR(nodeid) 0 0 GENTIMECURSOR(nodeid) 0 curdest 0 0];
        
        % Enqueue to the event list. In matlab you can append a row x to a
        % matrix X by using the command X=[X; x]
        elist1=[elist1; pkt];
    end

	% Return the CURTIME of node
	function t=getcurtime(node,pkt)
        if(pkt == 1)    
            idx=find(elist(:,SRC)==node,1,'first');
            t=elist(idx,CURTIME);
        else
            idx=find(elist1(:,SRC)==node,1,'first');
            t=elist1(idx,CURTIME);
        end        
        end
	
    function delaynodepkts(node,delay,pkt)
        if(pkt == 1)
              DELAYTIME=getcurtime(node, 1)+delay;
            % idx will have the row number of elist whose CURTIME field value
            % is less than DELAYTIME.
              list=find(elist(:,CURTIME)-DELAYTIME < 0 & elist(:,SRC)==node); 
            % Set the CURTIME field of all the rows in list to DELAYTIME.
              elist(list,CURTIME)=DELAYTIME;
        else
            DELAYTIME=getcurtime(node, 2)+delay;
            list1=find(elist1(:,CURTIME)-DELAYTIME < 0 & elist1(:,SRC)==node); 
            elist1(list1,CURTIME)=DELAYTIME;
        end    
	end
   
    % Move the first row of elist to SIMRESULT
    function updatesimlist()
        % elist(1,1:end) means row 1 and all columns.
        % In matlab you can append a row x to a
        % matrix X by using the command X=[X; x]
        % Enqueue to the the first row of event list to SIMRESULT. 
        SIMRESULT=[SIMRESULT; elist(1,1:end)];
        % This command will delete the first row of elist.
        elist(1,:)=[];
    end

    function updatesimlist1()
        % elist(1,1:end) means row 1 and all columns.
        % In matlab you can append a row x to a
        % matrix X by using the command X=[X; x]
        % Enqueue to the the first row of event list to SIMRESULT. 
        SIMRESULT1=[SIMRESULT1; elist1(1,1:end)];
        % This command will delete the first row of elist.
        elist1(1,:)=[];
    end

    function routingdelay = getRoutingPath(source, flag)
      if flag == 1
          initializeMatrix();
      end
      cost = 1000;
      router = 0;
      Row = M(source,:);
      for row = 1 : 4
          if Row(row) ~= 0 && Row(row) < cost
              cost = Row(row);
              router = row;
          end
      end
      elist(1,ROUTINGPATH) = router;
      routingdelay = tdelayr * cost;
      elist(1,PROPDELAY) = (2 * tdelay) + routingdelay;
    end

    function routingdelay = getRoutingPath1(source, flag)
      if flag == 1
          initializeMatrix();
      end
      cost = 1000;
      router = 0;
      Row = M(source,:);
      for row = 1 : 4
          if Row(row) ~= 0 && Row(row) < cost
              cost = Row(row);
              router = row;
          end
      end
      elist1(1,ROUTINGPATH) = router;
      routingdelay = tdelayr * cost;
      elist1(1,PROPDELAY) = (2 * tdelay) + routingdelay;
    end

    function initializeMatrix() 
      M(1,3) = round(unifrnd(1,10));
      M(3,1) = M(1,3);
      M(1,4) = round(unifrnd(1,10));
      M(4,1) = M(1,4);
      M(2,3) = round(unifrnd(1,10));
      M(3,2) = M(2,3);
      M(2,4) = round(unifrnd(1,10));  
      M(4,2) = M(2,4); 
      
      A3 = M(1,3) + M(3,2);
      A4 = M(1,4) + M(4,2);
      if A3 < A4
        M(1,2) = A3;
      else
        M(1,2) = A4;
      end
      
      B3 = M(2,3) + M(3,1);
      B4 = M(2,4) + M(4,1);
      if B3 < B4
        M(2,1) = B3;
      else
        M(2,1) = B4;
      end
      
      C1 = M(3,1) + M(1,4);
      C2 = M(3,2) + M(2,4);
      if C1 < C2
        M(3,4) = C1;
      else
        M(3,4) = C2;
      end
      
      D1 = M(4,1) + M(1,3);
      D2 = M(4,2) + M(2,3);
      if D1 < D2
        M(4,3) = D1;
      else
        M(4,3) = D2;
      end
    end

    function plotRoutingCounts()
        ArouteC = [SIMRESULT(:,SRC) == 1, SIMRESULT(:,ROUTINGPATH) == 3]; 
        ArouteCNum = length(ArouteC);
        ArouteCCnt = 0;
        for i = 1 : ArouteCNum
            if ArouteC(i,1) == 1 && ArouteC(i,2) == 1
                ArouteCCnt = ArouteCCnt + 1;
            end
        end
        
        ArouteD = [SIMRESULT(:,SRC) == 1, SIMRESULT(:,ROUTINGPATH) == 4]; 
        ArouteDNum = length(ArouteD);
        ArouteDCnt = 0;
        for i = 1 : ArouteDNum
            if ArouteD(i,1) == 1 && ArouteD(i,2) == 1
                ArouteDCnt = ArouteDCnt + 1;
            end
        end
    
        
        figure;
        subplot(2,2,1)
        X = [0,1,2,3,4,5];
        Y = [0,0,0,ArouteCCnt,ArouteDCnt,0];
        stem(X,Y);
        xlabel('Router Number');
        ylabel('Number of packets');
        title('Routing path taken by packets from A');

        BrouteC = [SIMRESULT(:,SRC) == 2, SIMRESULT(:,ROUTINGPATH) == 3]; 
        BrouteCNum = length(BrouteC);
        BrouteCCnt = 0;
        for i = 1 : BrouteCNum
            if BrouteC(i,1) == 1 && BrouteC(i,2) == 1
                BrouteCCnt = BrouteCCnt + 1;
            end
        end
        BrouteD = [SIMRESULT(:,SRC) == 2, SIMRESULT(:,ROUTINGPATH) == 4]; 
        BrouteDNum = length(BrouteD);
        BrouteDCnt = 0;
        for i = 1 : BrouteDNum
            if BrouteD(i,1) == 1 && BrouteD(i,2) == 1
                BrouteDCnt = BrouteDCnt + 1;
            end
        end
        
        subplot(2,2,2)
        X = [0,1,2,3,4,5];
        Y = [0,0,0,BrouteCCnt,BrouteDCnt,0];
        stem(X,Y);
        xlabel('Router Number');
        ylabel('Number of packets');
        title('Routing path taken by packets from B');

        CrouteA = [SIMRESULT1(:,SRC) == 3, SIMRESULT1(:,ROUTINGPATH) == 1]; 
        CrouteANum = length(CrouteA);
        CrouteACnt = 0;
        for i = 1 : CrouteANum
            if CrouteA(i,1) == 1 && CrouteA(i,2) == 1
                CrouteACnt = CrouteACnt + 1;
            end
        end      
        CrouteB = [SIMRESULT1(:,SRC) == 3, SIMRESULT1(:,ROUTINGPATH) == 2]; 
        CrouteBNum = length(CrouteB);
        CrouteBCnt = 0;
        for i = 1 : CrouteBNum
            if CrouteB(i,1) == 1 && CrouteB(i,2) == 1
                CrouteBCnt = CrouteBCnt + 1;
            end
        end      
        
        subplot(2,2,3)
        X = [0,1,2,3,4,5];
        Y = [0,CrouteACnt,CrouteBCnt,0,0,0];
        stem(X,Y);
        xlabel('Router Number');
        ylabel('Number of packets');
        title('Routing path taken by packets from C');

        DrouteA = [SIMRESULT1(:,SRC) == 4, SIMRESULT1(:,ROUTINGPATH) == 1]; 
        DrouteANum = length(DrouteA);
        DrouteACnt = 0;
        for i = 1 : DrouteANum
            if DrouteA(i,1) == 1 && DrouteA(i,2) == 1
                DrouteACnt = DrouteACnt + 1;
            end
        end
        DrouteB = [SIMRESULT1(:,SRC) == 4, SIMRESULT1(:,ROUTINGPATH) == 2]; 
        DrouteBNum = length(DrouteB);
        DrouteBCnt = 0;
        for i = 1 : DrouteBNum
            if DrouteB(i,1) == 1 && DrouteB(i,2) == 1
                DrouteBCnt = DrouteBCnt + 1;
            end
        end
        
        subplot(2,2,4)
        X = [0,1,2,3,4,5];
        Y = [0,DrouteACnt,DrouteBCnt,0,0,0];
        stem(X,Y);
        xlabel('Router Number');
        ylabel('Number of packets');
        title('Routing path taken by packets from D');
    end

    function plotPropDelayFromAandB()
        ACList = [];
        AproptoC = [SIMRESULT(:,SRC) == 1, SIMRESULT(:,DESTINATION) == 3, SIMRESULT(:,PROPDELAY)];
        ApropCNum = length(AproptoC);
        ApropCCnt = 0;
        for i = 1 : ApropCNum
            if AproptoC(i,1) == 1 && AproptoC(i,2) == 1
                propdelayAC = AproptoC(i,3);
                ApropCCnt = ApropCCnt + 1;
                rowAC = [ApropCCnt propdelayAC];
                ACList = [ACList; rowAC];
            end
        end
        
        figure;
        subplot(2,2,1)
        ACListpd = ACList(:,2);
        plot(1:ApropCCnt, ACListpd(1:end));
        xlabel('Packet Sequence #');
        ylabel('Delay in \mu sec');
        title('Delay from node A to C');
 
        ADList = [];
        AproptoD = [SIMRESULT(:,SRC) == 1, SIMRESULT(:,DESTINATION) == 4, SIMRESULT(:,PROPDELAY)];
        ApropDNum = length(AproptoD);
        ApropDCnt = 0;
        for i = 1 : ApropDNum
            if AproptoD(i,1) == 1 && AproptoD(i,2) == 1
                propdelayAD = AproptoD(i,3);
                ApropDCnt = ApropDCnt + 1;
                rowAD = [ApropDCnt propdelayAD];
                ADList = [ADList; rowAD];
            end
        end
        
        subplot(2,2,2)
        ADListpd = ADList(:,2);
        plot(1:ApropDCnt, ADListpd(1:end));
        xlabel('Packet Sequence #');
        ylabel('Delay in \mu sec');
        title('Delay from node A to D');

        BCList = [];
        BproptoC = [SIMRESULT(:,SRC) == 2, SIMRESULT(:,DESTINATION) == 3, SIMRESULT(:,PROPDELAY)];
        BpropCNum = length(BproptoC);
        BpropCCnt = 0;
        for i = 1 : BpropCNum
            if BproptoC(i,1) == 1 && BproptoC(i,2) == 1
                propdelayBC = BproptoC(i,3);
                BpropCCnt = BpropCCnt + 1;
                rowBC = [BpropCCnt propdelayBC];
                BCList = [BCList; rowBC];
            end
        end
        
        subplot(2,2,3)
        BCListpd = BCList(:,2);
        plot(1:BpropCCnt, BCListpd(1:end));
        xlabel('Packet Sequence #');
        ylabel('Delay in \mu sec');
        title('Delay from node B to C');

        BDList = [];
        BproptoD = [SIMRESULT(:,SRC) == 2, SIMRESULT(:,DESTINATION) == 4, SIMRESULT(:,PROPDELAY)];
        BpropDNum = length(BproptoD);
        BpropDCnt = 0;
        for i = 1 : BpropDNum
            if BproptoD(i,1) == 1 && BproptoD(i,2) == 1
                propdelayBD = BproptoD(i,3);
                BpropDCnt = BpropDCnt + 1;
                rowBD = [BpropDCnt propdelayBD];
                BDList = [BDList; rowBD];
            end
        end
        
        subplot(2,2,4)
        BDListpd = BDList(:,2);
        plot(1:BpropDCnt, BDListpd(1:end));
        xlabel('Packet Sequence #');
        ylabel('Delay in \mu sec');
        title('Delay from node B to D');
    end

    function plotPropDelayFromCandD()
        CAList = [];
        CproptoA = [SIMRESULT1(:,SRC) == 3, SIMRESULT1(:,DESTINATION) == 1, SIMRESULT1(:,PROPDELAY)];
        CpropANum = length(CproptoA);
        CpropACnt = 0;
        for i = 1 : CpropANum
            if CproptoA(i,1) == 1 && CproptoA(i,2) == 1
                propdelayCA = CproptoA(i,3);
                CpropACnt = CpropACnt + 1;
                rowCA = [CpropACnt propdelayCA];
                CAList = [CAList; rowCA];
            end
        end
        
        figure;
        subplot(2,2,1)
        CAListpd = CAList(:,2);
        plot(1:CpropACnt, CAListpd(1:end));
        xlabel('Packet Sequence #');
        ylabel('Delay in \mu sec');
        title('Delay from node C to A');
 
        CBList = [];
        CproptoB = [SIMRESULT1(:,SRC) == 3, SIMRESULT1(:,DESTINATION) == 2, SIMRESULT1(:,PROPDELAY)];
        CpropBNum = length(CproptoB);
        CpropBCnt = 0;
        for i = 1 : CpropBNum
            if CproptoB(i,1) == 1 && CproptoB(i,2) == 1
                propdelayCB = CproptoB(i,3);
                CpropBCnt = CpropBCnt + 1;
                rowCB = [CpropBCnt propdelayCB];
                CBList = [CBList; rowCB];
            end
        end
        
        subplot(2,2,2)
        CBListpd = CBList(:,2);
        plot(1:CpropBCnt, CBListpd(1:end));
        xlabel('Packet Sequence #');
        ylabel('Delay in \mu sec');
        title('Delay from node C to B');

        DAList = [];
        DproptoA = [SIMRESULT1(:,SRC) == 4, SIMRESULT1(:,DESTINATION) == 1, SIMRESULT1(:,PROPDELAY)];
        DpropANum = length(DproptoA);
        DpropACnt = 0;
        for i = 1 : DpropANum
            if DproptoA(i,1) == 1 && DproptoA(i,2) == 1
                propdelayDA = DproptoA(i,3);
                DpropACnt = DpropACnt + 1;
                rowDA = [DpropACnt propdelayDA];
                DAList = [DAList; rowDA];
            end
        end
        
        subplot(2,2,3)
        DAListpd = DAList(:,2);
        plot(1:DpropACnt, DAListpd(1:end));
        xlabel('Packet Sequence #');
        ylabel('Delay in \mu sec');
        title('Delay from node D to A');

        DBList = [];
        DproptoB = [SIMRESULT1(:,SRC) == 4, SIMRESULT1(:,DESTINATION) == 2, SIMRESULT1(:,PROPDELAY)];
        DpropBNum = length(DproptoB);
        DpropBCnt = 0;
        for i = 1 : DpropBNum
            if DproptoB(i,1) == 1 && DproptoB(i,2) == 1
                propdelayDB = DproptoB(i,3);
                DpropBCnt = DpropBCnt + 1;
                rowDB = [DpropBCnt propdelayDB];
                DBList = [DBList; rowDB];
            end
        end
        
        subplot(2,2,4)
        DBListpd = DBList(:,2);
        plot(1:DpropBCnt, DBListpd(1:end));
        xlabel('Packet Sequence #');
        ylabel('Delay in \mu sec');
        title('Delay from node D to B');
    end




disp(toc);

end
