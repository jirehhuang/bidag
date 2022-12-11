# bcb version of scorepossibleparents.PLUS1()

bcb_scorepossibleparents.PLUS1 <- function(parenttable,
                                           plus1lists,
                                           n,
                                           param,
                                           updatenodes,
                                           parentmaps,
                                           numparents,
                                           numberofparentsvec){

  listy <- vector("list", n)
  aliases <- plus1lists$aliases

  ## for every node which needs to be updated
  ## computing scores for node i?
  for (i in updatenodes){

    k <- nrow(aliases[[i]])
    ncols <- ncol(aliases[[i]])
    listz <- vector("list", k)

    ## compute score from observations for which i is not intervened
    temp_data <- param$data[param$settings$interventions != param$labels[i],,
                            drop = FALSE]
    temp_param <- BiDAG::scoreparameters(scoretype = "bde",
                                         data = temp_data)

    ## for every list
    ## computing scores for adding j -> i?
    for (j in seq_len(k)){

      if (j == 1){

        scoretemp <-
          BiDAG:::TableDAGscore.alias(parentrows = parenttable[[i]],
                                      j = i,
                                      n = n,
                                      alias = aliases[[i]][j, which(!is.na(aliases[[i]][j,]))],
                                      param = temp_param,
                                      parentmaps = parentmaps[[i]],
                                      numparents = numparents[i],
                                      numberofparentsvec = numberofparentsvec[[i]])
      } else{

        scoretemp <-
          BiDAG:::TableDAGscore.alias.plus1(parentrows = parenttable[[i]],
                                            j = i,
                                            n = n,
                                            alias = aliases[[i]][j,],
                                            param = temp_param,
                                            parentmaps = parentmaps[[i]],
                                            numparents = numparents[i],
                                            numberofparentsvec = numberofparentsvec[[i]])
      }
      listz[[j]] <- as.matrix(scoretemp)
    }
    listy[[i]] <- listz
  }
  return(listy)
}



## bcb version of TableDAGscore.alias() and TableDAGscore.alias.plus1()

bcb_TableDAGscore.alias <- function(parentrows,
                                    j,
                                    n,
                                    alias,
                                    param,
                                    parentmaps = NULL,
                                    numparents = NULL,
                                    numberofparentsvec = NULL,
                                    plus1 = FALSE){

  nrows <- nrow(parentrows)
  if (plus1){

    parentnodes <- alias[parentrows[nrows, !is.na(parentrows[nrows,])] + 1]
    addpar <- alias[1]

  } else{

    parentnodes <- alias[parentrows[nrows, !is.na(parentrows[nrows,])]]
  }
  debug <- 0  # TODO: add debugging
  network <- bnlearn::empty.graph(nodes = param$labels)
  amat <- bnlearn::amat(network)

  scrs <- apply(parentrows, 1, function(x){

    # direct parents -> j in network
    if (plus1){

      amat[union(parentnodes[x], addpar), j] <- 1

    } else{

      amat[parentnodes[x], j] <- 1
    }
    bnlearn::amat(network) <- amat

    ## TODO: Gaussian case
    extra.args <- list(iss = param$chi)
    suppressWarnings(
      extra.args <- bnlearn:::check.score.args(score = param$settings$score, network = network,
                                               data = param$settings$data, extra.args = extra.args)
    )
    ## compute and store score
    scr <- local_score(network = network, data = param$settings$data, score = param$settings$score,
                       targets = param$settings$nodes[j], extra.args = extra.args,
                       interventions = param$settings$interventions, debug = debug >= 4)
    return(scr)
  })
  return(scrs)
}



# bcb version of TableDAGscore.alias.plus1()

bcb_TableDAGscore.alias.plus1 <- function(parentrows,
                                          j,
                                          n,
                                          alias,
                                          param,
                                          parentmaps = NULL,
                                          numparents = NULL,
                                          numberofparentsvec = NULL){

  nrows <- nrow(parentrows)
  parentnodes <- alias[parentrows[nrows, !is.na(parentrows[nrows,])] + 1]
  addpar <- alias[1]

  debug <- 0  # TODO: add debugging
  network <- bnlearn::empty.graph(nodes = param$labels)
  amat <- bnlearn::amat(network)

  scrs <- apply(parentrows, 1, function(x){

    # direct parents -> j in network
    amat[union(parentnodes[x], addpar), j] <- 1
    bnlearn::amat(network) <- amat

    ## TODO: Gaussian case
    extra.args <- list(iss = param$chi)
    suppressWarnings(
      extra.args <- bnlearn:::check.score.args(score = param$settings$score, network = network,
                                               data = param$data, extra.args = extra.args)
    )
    ## compute and store score
    scr <- local_score(network = network, data = param$data0, score = param$settings$score,
                       targets = param$settings$nodes[j], extra.args = extra.args,
                       interventions = param$interventions, debug = debug >= 4)
    return(scr)
  })
  return(scrs)
}
