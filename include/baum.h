#ifdef __cplusplus
extern "C" {
#endif



#ifndef _BAUM_H
#define _BAUM_H

#include "list.h"

#define TREE_HEAD(name)						\
struct tree_node name = {&name,					\
			 { NULL, NULL},				\
			 {&name.children, &name.children}	\
}

#define INIT_TREE_HEAD(ptr) do {		\
	(ptr)->parent = (ptr);			\
	INIT_LIST_HEAD(&(ptr)->children);	\
} while(0);

struct tree_node{
	struct tree_node *parent;
  	struct list_head node;
	struct list_head children;
};

#define tree_entry(ptr, type, member) \
	((type *)((char *)(ptr)-(unsigned long)(&((type *)0)->member)))

#define child_entry(ptr, type, member) \
	tree_entry(list_entry(ptr, struct tree_node, node), type, member)

#define first_child(root, member) \
	(root)->member.children.next

#if 0
struct tree_node {
	struct tree_node *parent;
	struct tree_node *prev, *next;
	struct tree_node *head_prev, *head_next;
	struct my_funny_data_struct data;
}
#endif

void delete_tree_node(struct tree_node *old);
void insert_tree_node(struct tree_node *parent, struct tree_node *new);

#endif


#ifdef __cplusplus
}
#endif

