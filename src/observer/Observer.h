/*  PartitionTest, fast selection of the best fit partitioning scheme for
 *  multi-gene data sets.
 *  Copyright May 2013 by Diego Darriba
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other inquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file Observer.h
 */

#ifndef OBSERVER_H_
#define OBSERVER_H_

#include "model/Model.h"
#include "model/ModelSet.h"
#include "options/ParTestOptions.h"
#include <time.h>
#include <assert.h>

namespace partest {

/**
 * @brief Type of the message sent by the observable instances.
 */
enum MessageType {
  MT_GROUP_INIT, /** Start of a group of schemes computation */
  MT_SCHEME_INIT, /** Start of a partitioning scheme computation */
  MT_MODELSET_INIT, /** Start of a model set computation */
  MT_FTREE_INIT,
  MT_FTREE_END,
  MT_SINGLE_INIT,
  MT_SINGLE_END,
  MT_MODELSET_END, /** End of a model set computation */
  MT_SCHEME_END, /** End of a partitioning scheme computation */
  MT_GROUP_END /** End of a group of schemes computation */
};

/**
 * @brief Structure for communications between partest::Observer and partest::Observable instances.
 */
struct ObservableInfo {
	/**
	 * @brief Instantiates a new partest::ObservableInfo.
	 *
	 * @param[in] type Message type.
	 * @param[in] p_index Index with some meaning for the observers.
	 * @param[in] model Reference to a model (if appropriate).
	 * @param[in] time Timestamp with some meaning for the observers.
	 * @param[in] current_index Current index within the group.
	 * @param[in] max_index Maximum index within the group.
	 * @param[in] message Human readable message.
	 */
    ObservableInfo(MessageType type, t_partitionElementId p_index, Model * model, time_t time = 0,
        t_partitionElementId current_index = 0, t_partitionElementId max_index = 0, string message = 0) :
        type(type), model(model), time(time), current_index(current_index), max_index(
            max_index), p_index(p_index), message(message) {
    	assert(type != MT_MODELSET_INIT && type != MT_MODELSET_END);
    	modelset = 0;
    }
    /**
    	 * @brief Instantiates a new partest::ObservableInfo.
    	 *
    	 * @param[in] type Message type.
    	 * @param[in] p_index Index with some meaning for the observers.
    	 * @param[in] modelset Reference to a model set (if appropriate).
    	 * @param[in] time Timestamp with some meaning for the observers.
    	 * @param[in] current_index Current index within the group.
    	 * @param[in] max_index Maximum index within the group.
    	 * @param[in] message Human readable message.
    	 */
    ObservableInfo(MessageType type, t_partitionElementId p_index, ModelSet * modelset, time_t time = 0,
    		t_partitionElementId current_index = 0, t_partitionElementId max_index = 0, string message = 0) :
        type(type), modelset(modelset), time(time), current_index(current_index), max_index(
            max_index), p_index(p_index), message(message) {
    	assert(type == MT_MODELSET_INIT || type == MT_MODELSET_END);
    	model = 0;
    }
    const MessageType type; /** Type of the message sent by the observable instances. */
    Model * model; /** Reference to a model (if appropriate). */
    ModelSet * modelset; /** Reference to a model set (if appropriate). */
    const time_t time; /** Timestamp with some meaning for the observers. */
    const t_partitionElementId current_index; /** Current index within the group. */
    const t_partitionElementId max_index; /** Maximum index within the group. */
    const t_partitionElementId p_index; /** Index with some meaning for the observers. */
    string message; /** Human readable message. */
};

/**
 * @brief Interface for observer classes.
 */
class Observer {
  public:
    Observer() {}
    virtual ~Observer() {}
    /**
     * @brief Message for get a notification from the partest::Observable instances.
     *
     * @param info Message from the partest::Observable
     * @param run_instance Optimization options.
     */
    virtual void update(const ObservableInfo & info,
        ParTestOptions * run_instance = NULL)= 0;
};

} /* namespace partest */
#endif /* OBSERVER_H_ */
