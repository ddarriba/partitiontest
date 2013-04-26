/*
 * Observer.hpp
 *
 *  Created on: 30/05/2012
 *      Author: diego
 */

#ifndef OBSERVER_H_
#define OBSERVER_H_

#include "../model/Model.h"
#include "../model/ModelSet.h"
#include "../options/ParTestOptions.h"
#include <time.h>
#include <assert.h>

namespace partest {

enum MessageType {
  MT_GROUP_INIT,
  MT_PARTITION_INIT,
  MT_MODELSET_INIT,
  MT_FTREE_INIT,
  MT_FTREE_END,
  MT_SINGLE_INIT,
  MT_SINGLE_END,
  MT_MODELSET_END,
  MT_PARTITION_END,
  MT_GROUP_END
};
struct ObservableInfo {
    ObservableInfo(MessageType type, unsigned int p_index, Model * model, time_t time = 0,
        unsigned int current_index = 0, unsigned int max_index = 0, string * message = 0) :
        type(type), model(model), time(time), current_index(current_index), max_index(
            max_index), p_index(p_index), message(message) {
    	assert(type != MT_MODELSET_INIT && type != MT_MODELSET_END);
    	modelset = 0;
    }
    ObservableInfo(MessageType type, unsigned int p_index, ModelSet * modelset, time_t time = 0,
        unsigned int current_index = 0, unsigned int max_index = 0, string * message = 0) :
        type(type), modelset(modelset), time(time), current_index(current_index), max_index(
            max_index), p_index(p_index), message(message) {
    	assert(type == MT_MODELSET_INIT || type == MT_MODELSET_END);
    	model = 0;
    }
    const MessageType type;
    Model * model;
    ModelSet * modelset;
    const time_t time;
    const unsigned int current_index;
    const unsigned int max_index;
    const unsigned int p_index;
    string * message;
//    const t_partition partition;
};

class Observer {
  public:
    Observer() {}
    virtual ~Observer() {}
    virtual void update(const ObservableInfo & info,
        ParTestOptions * run_instance = NULL)= 0;
};

} /* namespace partest */
#endif /* OBSERVER_H_ */
